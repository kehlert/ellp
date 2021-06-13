use std::io::Read;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ellp::parse_mps;
use ellp::Problem;

fn fibonacci(n: u64) -> u64 {
    let f = std::fs::File::open("prob.mps").unwrap();
    let mut reader = std::io::BufReader::new(f);
    let mut mps = Vec::new();
    reader.read_to_end(&mut mps).unwrap();
    let prob: Problem = parse_mps(mps.as_slice()).unwrap();

    match n {
        0 => 1,
        1 => 1,
        n => fibonacci(n - 1) + fibonacci(n - 2),
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("fib 20", |b| b.iter(|| fibonacci(black_box(20))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
