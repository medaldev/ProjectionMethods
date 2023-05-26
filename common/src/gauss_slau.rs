pub fn gauss(A: &mut Vec<Vec<f64>>) -> Vec<f64>{

    let n: usize = A.len();

    for k in 0..n {
        for j in k+1..n {
            let dcoeff: f64 = A[j][k] / A[k][k];
            for i in k..n {
                A[j][i] -= dcoeff * A[k][i];
            }
            A[j][n] -= dcoeff * A[k][n];


        }
        for j in k+1..n {
            let dcoeff: f64 = A[n-j-1][n-k-1] / A[n-k-1][n-k-1];
            for i in k..n {
                A[n-j-1][n-i-1] -= dcoeff * A[n-k-1][n-i-1];
            }
            A[n-j-1][n] -= dcoeff * A[n-k-1][n];

        }
    }

    let mut c: Vec<f64> = Vec::with_capacity(n + 1);

    for i in 0..n {
        c.push(A[i][n] / A[i][i]);
    }

    c
}
