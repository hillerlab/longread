//! Integration tests for `longread rg` (normalize) and `longread check` (CCS-chunk validation).
//!
//! Tiny BAMs are built in-memory with `noodles`, written to a temp dir, processed through the
//! public `pacbio` API, and read back for assertions.

use std::fs;
use std::path::Path;

use noodles_bam as bam;
use noodles_sam as sam;
use sam::alignment::io::Write as _;
use sam::alignment::record::data::field::Tag;
use sam::alignment::record_buf::data::field::Value;
use sam::alignment::RecordBuf;
use sam::header::record::value::map::{read_group, ReadGroup};
use sam::header::record::value::Map;

use longread::pacbio::check::{self, CheckParams};
use longread::pacbio::normalize::{self, NormalizeParams};

const ZM: Tag = Tag::new(b'z', b'm');

/// Build a subread/CCS BAM header with a single read group.
fn header(rg_id: &str, readtype: &str, platform_unit: &str) -> sam::Header {
    let mut rg = Map::<ReadGroup>::default();
    rg.other_fields_mut().insert(
        read_group::tag::DESCRIPTION,
        format!("READTYPE={readtype};BINDINGKIT=100")
            .into_bytes()
            .into(),
    );
    rg.other_fields_mut().insert(
        read_group::tag::PLATFORM_UNIT,
        platform_unit.as_bytes().to_vec().into(),
    );
    sam::Header::builder().add_read_group(rg_id, rg).build()
}

/// Build one unmapped record with the given QNAME and optional `RG`/`zm` tags.
fn record(qname: &str, rg: Option<&str>, zm: Option<i32>) -> RecordBuf {
    let mut rec = RecordBuf::default();
    *rec.name_mut() = Some(qname.as_bytes().to_vec().into());
    if let Some(rg) = rg {
        rec.data_mut().insert(
            Tag::READ_GROUP,
            Value::String(rg.as_bytes().to_vec().into()),
        );
    }
    if let Some(zm) = zm {
        rec.data_mut().insert(ZM, Value::Int32(zm));
    }
    rec
}

fn write_bam(path: &Path, header: &sam::Header, records: &[RecordBuf]) {
    let mut writer = bam::io::Writer::new(fs::File::create(path).unwrap());
    writer.write_header(header).unwrap();
    for rec in records {
        writer.write_alignment_record(header, rec).unwrap();
    }
    writer.try_finish().unwrap();
}

fn read_bam(path: &Path) -> (sam::Header, Vec<RecordBuf>) {
    let mut reader = bam::io::Reader::new(fs::File::open(path).unwrap());
    let header = reader.read_header().unwrap();
    let mut records = Vec::new();
    let mut rec = RecordBuf::default();
    while reader.read_record_buf(&header, &mut rec).unwrap() != 0 {
        records.push(rec.clone());
    }
    (header, records)
}

fn name_of(rec: &RecordBuf) -> String {
    rec.name().unwrap().to_string()
}

fn int_tag(rec: &RecordBuf, tag: Tag) -> Option<i64> {
    match rec.data().get(&tag)? {
        Value::Int8(v) => Some(i64::from(*v)),
        Value::UInt8(v) => Some(i64::from(*v)),
        Value::Int16(v) => Some(i64::from(*v)),
        Value::UInt16(v) => Some(i64::from(*v)),
        Value::Int32(v) => Some(i64::from(*v)),
        Value::UInt32(v) => Some(i64::from(*v)),
        _ => None,
    }
}

fn string_tag(rec: &RecordBuf, tag: Tag) -> Option<String> {
    match rec.data().get(&tag)? {
        Value::String(s) => Some(s.to_string()),
        _ => None,
    }
}

/// A default subread run: two source movies with overlapping ZMW numbering.
fn write_two_movie_inputs(dir: &Path) -> (std::path::PathBuf, std::path::PathBuf) {
    let a = dir.join("movie.a.subreads.bam");
    let b = dir.join("movie.b.subreads.bam");
    // movie.a: ZMW 5,6,7 (two subreads on ZMW 5)
    write_bam(
        &a,
        &header("ffffffff", "SUBREAD", "movie.raw"),
        &[
            record("movie.a/5/0_10", Some("ffffffff"), Some(5)),
            record("movie.a/5/12_22", Some("ffffffff"), Some(5)),
            record("movie.a/6/0_10", Some("ffffffff"), Some(6)),
            record("movie.a/7/0_10", Some("ffffffff"), Some(7)),
        ],
    );
    // movie.b: ZMW 1,2
    write_bam(
        &b,
        &header("ffffffff", "SUBREAD", "movie.raw"),
        &[
            record("movie.b/1/0_10", Some("ffffffff"), Some(1)),
            record("movie.b/2/0_10", Some("ffffffff"), Some(2)),
        ],
    );
    (a, b)
}

fn normalize_params(bams: Vec<std::path::PathBuf>, dir: &Path, threads: usize) -> NormalizeParams {
    NormalizeParams {
        bams,
        movie: "movie.test".to_string(),
        outdir: dir.to_path_buf(),
        zmw_map: dir.join("zmw_map.tsv"),
        threads,
        max_open_files: 64,
    }
}

#[test]
fn normalize_assigns_disjoint_global_zmws() {
    let dir = tempfile::tempdir().unwrap();
    let (a, b) = write_two_movie_inputs(dir.path());

    let stats = normalize::run(&normalize_params(vec![a, b], dir.path(), 0)).unwrap();
    assert_eq!(stats.input_bams, 2);
    assert_eq!(stats.source_files, 2);
    assert_eq!(stats.records, 6);
    // movie.a span [5,7] -> 3, movie.b span [1,2] -> 2; total capacity 5.
    assert_eq!(stats.zmw_capacity, 5);
    assert_eq!(stats.outputs.len(), 2);

    // movie.a sorts before movie.b: a offset 0, b offset 3.
    // a: 5->1, 6->2, 7->3 ; b: 1->4, 2->5.
    let (ha, ra) = read_bam(&dir.path().join("movie.a.subreads.normalized.bam"));
    assert_eq!(ha.read_groups().len(), 1);
    let (key, map) = ha.read_groups().first().unwrap();
    assert_eq!(key.to_string(), stats.rg_id);
    assert_eq!(
        map.other_fields()
            .get(&read_group::tag::PLATFORM_UNIT)
            .unwrap()
            .to_string(),
        "movie.test"
    );

    let names: Vec<String> = ra.iter().map(name_of).collect();
    assert_eq!(
        names,
        vec![
            "movie.test/1/0_10",
            "movie.test/1/12_22",
            "movie.test/2/0_10",
            "movie.test/3/0_10",
        ]
    );
    for rec in &ra {
        assert_eq!(
            string_tag(rec, Tag::READ_GROUP).as_deref(),
            Some(stats.rg_id.as_str())
        );
    }
    assert_eq!(int_tag(&ra[0], ZM), Some(1));
    assert_eq!(int_tag(&ra[2], ZM), Some(2));
    assert_eq!(int_tag(&ra[3], ZM), Some(3));

    let (_, rb) = read_bam(&dir.path().join("movie.b.subreads.normalized.bam"));
    let names_b: Vec<String> = rb.iter().map(name_of).collect();
    assert_eq!(names_b, vec!["movie.test/4/0_10", "movie.test/5/0_10"]);
    assert_eq!(int_tag(&rb[0], ZM), Some(4));
    assert_eq!(int_tag(&rb[1], ZM), Some(5));
}

#[test]
fn normalize_writes_expected_zmw_map() {
    let dir = tempfile::tempdir().unwrap();
    let (a, b) = write_two_movie_inputs(dir.path());
    normalize::run(&normalize_params(vec![a, b], dir.path(), 0)).unwrap();

    let map = fs::read_to_string(dir.path().join("zmw_map.tsv")).unwrap();
    let expected = "source_bam\toriginal_movie\toriginal_zmw\tmovie\tzmw\n\
                    movie.a.subreads.bam\tmovie.a\t5\tmovie.test\t1\n\
                    movie.a.subreads.bam\tmovie.a\t6\tmovie.test\t2\n\
                    movie.a.subreads.bam\tmovie.a\t7\tmovie.test\t3\n\
                    movie.b.subreads.bam\tmovie.b\t1\tmovie.test\t4\n\
                    movie.b.subreads.bam\tmovie.b\t2\tmovie.test\t5\n";
    assert_eq!(map, expected);
}

#[test]
fn normalize_is_byte_identical_across_thread_counts() {
    let dir1 = tempfile::tempdir().unwrap();
    let dir2 = tempfile::tempdir().unwrap();
    let (a1, b1) = write_two_movie_inputs(dir1.path());
    let (a2, b2) = write_two_movie_inputs(dir2.path());

    normalize::run(&normalize_params(vec![a1, b1], dir1.path(), 1)).unwrap();
    normalize::run(&normalize_params(vec![a2, b2], dir2.path(), 8)).unwrap();

    for name in [
        "movie.a.subreads.normalized.bam",
        "movie.b.subreads.normalized.bam",
        "zmw_map.tsv",
    ] {
        let x = fs::read(dir1.path().join(name)).unwrap();
        let y = fs::read(dir2.path().join(name)).unwrap();
        assert_eq!(x, y, "{name} differs between 1 and 8 threads");
    }
}

#[test]
fn normalize_rejects_non_subread_read_group() {
    let dir = tempfile::tempdir().unwrap();
    let bam = dir.path().join("movie.c.subreads.bam");
    write_bam(
        &bam,
        &header("ffffffff", "CCS", "movie.raw"),
        &[record("movie.c/1/ccs", Some("ffffffff"), Some(1))],
    );
    let err = normalize::run(&normalize_params(vec![bam], dir.path(), 0)).unwrap_err();
    assert!(err.to_string().contains("READTYPE=SUBREAD"), "{err}");
}

#[test]
fn normalize_rejects_missing_zm_tag() {
    let dir = tempfile::tempdir().unwrap();
    let bam = dir.path().join("movie.d.subreads.bam");
    write_bam(
        &bam,
        &header("ffffffff", "SUBREAD", "movie.raw"),
        &[record("movie.d/1/0_10", Some("ffffffff"), None)],
    );
    let err = normalize::run(&normalize_params(vec![bam], dir.path(), 0)).unwrap_err();
    assert!(err.to_string().contains("lacks RG or zm"), "{err}");
}

#[test]
fn normalize_rejects_missing_rg_tag() {
    let dir = tempfile::tempdir().unwrap();
    let bam = dir.path().join("movie.e.subreads.bam");
    write_bam(
        &bam,
        &header("ffffffff", "SUBREAD", "movie.raw"),
        &[record("movie.e/1/0_10", None, Some(1))],
    );
    let err = normalize::run(&normalize_params(vec![bam], dir.path(), 0)).unwrap_err();
    assert!(err.to_string().contains("lacks RG or zm"), "{err}");
}

#[test]
fn normalize_rejects_duplicate_qname_within_file() {
    let dir = tempfile::tempdir().unwrap();
    let bam = dir.path().join("chunk.subreads.bam");
    let h = header("ffffffff", "SUBREAD", "movie.raw");
    write_bam(
        &bam,
        &h,
        &[
            record("movie.y/1/0_10", Some("ffffffff"), Some(1)),
            record("movie.y/1/0_10", Some("ffffffff"), Some(1)),
        ],
    );
    let err = normalize::run(&normalize_params(vec![bam], dir.path(), 0)).unwrap_err();
    assert!(
        err.to_string().contains("duplicate QNAME within file"),
        "{err}"
    );
}

#[test]
fn normalize_disambiguates_colliding_movie_across_files() {
    // The `wgs` failure mode: two PBSIM3 files share a movie name and both restart ZMW at 1, so
    // the same original QNAME appears in both. Per-file allocation must give them distinct global
    // ZMWs rather than rejecting the input.
    let dir = tempfile::tempdir().unwrap();
    let f1 = dir.path().join("chunk1.subreads.bam");
    let f2 = dir.path().join("chunk2.subreads.bam");
    let h = header("ffffffff", "SUBREAD", "movie.raw");
    write_bam(
        &f1,
        &h,
        &[record("movie.x/1/0_10", Some("ffffffff"), Some(1))],
    );
    write_bam(
        &f2,
        &h,
        &[record("movie.x/1/0_10", Some("ffffffff"), Some(1))],
    );

    let stats = normalize::run(&normalize_params(vec![f1, f2], dir.path(), 0)).unwrap();
    assert_eq!(stats.source_files, 2);
    assert_eq!(stats.zmw_capacity, 2);

    // chunk1 sorts before chunk2 -> global ZMW 1 and 2 respectively.
    let (_, r1) = read_bam(&dir.path().join("chunk1.subreads.normalized.bam"));
    let (_, r2) = read_bam(&dir.path().join("chunk2.subreads.normalized.bam"));
    assert_eq!(name_of(&r1[0]), "movie.test/1/0_10");
    assert_eq!(name_of(&r2[0]), "movie.test/2/0_10");
    assert_eq!(int_tag(&r1[0], ZM), Some(1));
    assert_eq!(int_tag(&r2[0], ZM), Some(2));
}

#[test]
fn normalize_rejects_malformed_qname() {
    let dir = tempfile::tempdir().unwrap();
    let bam = dir.path().join("movie.f.subreads.bam");
    write_bam(
        &bam,
        &header("ffffffff", "SUBREAD", "movie.raw"),
        &[record("movie.f/1", Some("ffffffff"), Some(1))],
    );
    let err = normalize::run(&normalize_params(vec![bam], dir.path(), 0)).unwrap_err();
    assert!(
        err.to_string().contains("malformed PBSIM3 subread QNAME"),
        "{err}"
    );
}

fn check_params(bams: Vec<std::path::PathBuf>, out: &Path) -> CheckParams {
    CheckParams {
        bams,
        out: out.to_path_buf(),
        threads: 0,
        max_open_files: 64,
    }
}

#[test]
fn check_passes_for_disjoint_ccs_chunks() {
    let dir = tempfile::tempdir().unwrap();
    let c1 = dir.path().join("ccs1.bam");
    let c2 = dir.path().join("ccs2.bam");
    let h = header("rg0", "CCS", "movie.test");
    write_bam(
        &c1,
        &h,
        &[
            record("movie.test/1/ccs", Some("rg0"), Some(1)),
            record("movie.test/2/ccs", Some("rg0"), Some(2)),
        ],
    );
    write_bam(
        &c2,
        &h,
        &[
            record("movie.test/3/ccs", Some("rg0"), Some(3)),
            record("movie.test/4/ccs", Some("rg0"), Some(4)),
        ],
    );
    let out = dir.path().join("ccs_chunks.valid");
    let stats = check::run(&check_params(vec![c1, c2], &out)).unwrap();
    assert_eq!(stats.records, 4);
    assert_eq!(stats.distinct_keys, 4);
    assert!(out.exists());
}

#[test]
fn check_fails_for_overlapping_ccs_chunks() {
    let dir = tempfile::tempdir().unwrap();
    let c1 = dir.path().join("ccs1.bam");
    let c2 = dir.path().join("ccs2.bam");
    let h = header("rg0", "CCS", "movie.test");
    write_bam(&c1, &h, &[record("movie.test/2/ccs", Some("rg0"), Some(2))]);
    // movie.test/2 also appears in the second chunk -> overlap.
    write_bam(&c2, &h, &[record("movie.test/2/ccs", Some("rg0"), Some(2))]);
    let out = dir.path().join("ccs_chunks.valid");
    let err = check::run(&check_params(vec![c1, c2], &out)).unwrap_err();
    assert!(err.to_string().contains("CCS chunks overlap"), "{err}");
    assert!(err.to_string().contains("movie.test/2"), "{err}");
    assert!(!out.exists());
}
