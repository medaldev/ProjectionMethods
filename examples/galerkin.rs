use common::{gauss, get_step, grid, print_matrix};


fn main() {
    let n: usize = 20;
    let (a, b) = (0., 1.);
    let lam = 0.1;
    let nq: usize = 100;

    let C = galerkin(n, a, b, lam, nq);


    println!("{:?}", C);


}

fn K(x: f64, s: f64) -> f64 {
    x - s
}

fn phi(i: usize, j: usize) -> f64{
    return ((i == j) as i32) as f64;
}


fn galerkin(n: usize, a:f64, b:f64, lam: f64, nq: usize) -> Vec<f64>{
    let h = get_step(a, b, n);

    let X = grid(a, b, n);
    let mut A: Vec<Vec<f64>> = vec![vec![0.0; n + 1]; n];

    for i in 0..n {
        for j in 0..n {
            A[i][j] = phi(i, j) * h - lam * intg2((X[i], X[i+1]), (X[j], X[j+1]), nq);
        }
        A[i][n] = intg1((X[i], X[i+1]), nq, lam);
    }

    let C = gauss(&mut A.clone());
    C
}

fn u0(x: f64, lam: f64) -> f64{
    1.0 - lam * (x - 1.0 / 2.0)
}


fn intg1(AB: (f64, f64), nq: usize, lam: f64) -> f64{
    let mut sum: f64 = 0.;
    let h1 = get_step(AB.0, AB.1, nq);

    for i1 in 0..nq {
        sum +=u0(AB.0 + (i1 as f64 + 0.5) * h1, lam)
    }

    sum * h1

}

fn intg2(AB: (f64, f64), CD: (f64, f64), nq: usize) -> f64{
    let mut sum: f64 = 0.;
    let h1 = get_step(AB.0, AB.1, nq);
    let h2 = get_step(CD.0, CD.1, nq);

    for i1 in 0..nq {
        for i2 in 0.. nq {
            sum +=K(AB.0 + (i1 as f64 + 0.5) * h1, CD.0 + (i2 as f64 + 0.5) * h2)
        }
    }

    sum * h1 * h2

}




