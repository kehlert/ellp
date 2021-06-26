use std::io::Read;

use criterion::{criterion_group, criterion_main, Criterion};
use ellp::{parse_mps, MpsParsingError, Problem};

const BENCH_PROB_REL_DIR: &str = "benches/benchmark_problems";

#[derive(Debug)]
struct ProblemPaths {
    mps: std::path::PathBuf,
    solution: std::path::PathBuf,
}

impl ProblemPaths {
    fn new(prob_name: &str) -> Self {
        let mut prob_dir = std::env::current_dir().unwrap();
        prob_dir.push(BENCH_PROB_REL_DIR);
        prob_dir.push(prob_name);

        let mps = prob_dir.join(format!("{}.mps", prob_name));
        let solution = prob_dir.join(format!("{}_solution.csv", prob_name));

        Self { mps, solution }
    }
}

fn read_mps_file(p: &std::path::Path) -> Result<Problem, MpsParsingError> {
    let f = std::fs::File::open(p).unwrap();
    let mut reader = std::io::BufReader::new(f);
    let mut mps = Vec::new();
    reader.read_to_end(&mut mps).unwrap();
    let mps = std::str::from_utf8(&mps).unwrap();
    parse_mps(mps)
}

fn neos_5251015() {
    let prob_paths = ProblemPaths::new("neos-5251015_lp");

    println!("\nPATHS\n{:?}", prob_paths);

    let prob = read_mps_file(&prob_paths.mps).unwrap();
    println!("{}", prob);
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("fib 20", |b| b.iter(neos_5251015));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
