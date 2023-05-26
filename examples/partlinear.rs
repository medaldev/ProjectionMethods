use common::{gauss, get_step, grid, print_matrix};


fn main() {
    let n: usize = 2;
    let (a, b) = (0., 1.);
    let lam = 0.1;
    let nq: usize = 100;

    let C = partlinear(n, a, b, lam, nq);


    println!("{:?}", C);
}

fn K(s: f64, t: f64) -> f64 {
    s - t
}


fn partlinear(n: usize, a: f64, b: f64, lam: f64, nq: usize) -> Vec<f64> {
    let h = get_step(a, b, n);

    let X = grid(a, b, n);
    let mut A: Vec<Vec<f64>> = vec![vec![0.0; n + 1]; n];

    let super_phi = |x: f64, AB: (f64, f64)| -> f64 {
        if AB == (X[0], X[1]) {
            phi_0(x, AB)
        } else if AB == (X[n - 1], X[n]) {
            phi_n(x, AB)
        } else {
            phi_i(x, AB)
        }
    };

    for i in 0..n {
        for j in 0..n {
            A[i][j] = intg2(lam, h, (X[i], X[i + 1]), (X[j], X[j + 1]), nq, &super_phi);
        }
        A[i][n] = intg1(lam, (X[i], X[i + 1]), nq, &super_phi);
    }

    print_matrix(&A);

    let C = gauss(&mut A.clone());
    C
}

fn u0(x: f64, lam: f64) -> f64 {
    1.0 - lam * (x - 1.0 / 2.0)
}


fn phi_i(x: f64, AB: (f64, f64)) -> f64 {
    let c: f64 = AB.0 + (AB.1 - AB.0) / 2.;

    if x >= AB.0 && x <= c {
        (x - AB.0) / (c - AB.0)
    } else if x >= c && x <= AB.1 {
        (AB.1 - x) / (AB.1 - c)
    } else {
        0.
    }
}


fn phi_0(x: f64, AB: (f64, f64)) -> f64 {
    if x >= AB.0 && x <= AB.1 { (AB.1 - x) / (AB.1 - AB.0) } else { 0. }
}

fn phi_n(x: f64, AB: (f64, f64)) -> f64 {
    if x >= AB.0 && x <= AB.1 { (x - AB.0) / (AB.1 - AB.0) } else { 0. }
}

fn intg2(lam: f64, h: f64, AB: (f64, f64), CD: (f64, f64), nq: usize,
         phi: &dyn Fn(f64, (f64, f64)) -> f64) -> f64 {

    let mut sum: f64 = 0.;
    let h_AB = get_step(AB.0, AB.1, nq);
    let h_CD = get_step(CD.0, CD.1, nq);

    for i1 in 0..nq {
        let x_AB = AB.0 + (i1 as f64 + 0.5) * h_AB;

        let left = phi(x_AB, AB) * phi(x_AB, CD); // Questions

        for i2 in 0..nq {
            let x_CD = CD.0 + (i2 as f64 + 0.5) * h_CD;

            sum += -lam * K(x_AB, x_CD) * phi(x_AB, AB) * phi(x_CD, CD)
        }

        sum += left
    }

    sum * h_AB * h_CD
}

fn intg1(lam: f64, AB: (f64, f64), nq: usize, phi: &dyn Fn(f64, (f64, f64)) -> f64) -> f64 {
    let mut sum: f64 = 0.;
    let h1 = get_step(AB.0, AB.1, nq);

    for i1 in 0..nq {
        let x_AB = AB.0 + (i1 as f64 + 0.5) * h1;

        sum += u0(x_AB, lam) * phi(x_AB, AB)
    }

    sum * h1
}




