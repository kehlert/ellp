use std::io::Read;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ellp::Problem;
use std::convert::TryInto;

fn fibonacci(n: u64) -> u64 {
    let f = std::fs::File::open("prob.mps").unwrap();
    let mut reader = std::io::BufReader::new(f);
    let mut mps = Vec::new();
    reader.read_to_end(&mut mps).unwrap();
    let prob: Problem = mps.as_slice().try_into().unwrap();

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
