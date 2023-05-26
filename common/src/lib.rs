mod gauss_slau;

use std::fmt::Display;
pub use gauss_slau::gauss;


pub fn get_step(a: f64, b: f64, n: usize) -> f64{
    (b - a) / n as f64
}

pub fn grid(a: f64, b: f64, n: usize) -> Vec<f64>{
    let h = get_step(a, b, n);
    let mut res: Vec<f64> = vec!();
    for i in 0..n + 1 {
        res.push(a + i as f64 * h);
    }
    res
}

pub fn print_matrix<T>(matrix: &Vec<Vec<T>>) where T: Display {
    for i in 0..matrix.len() {
        for j in 0..matrix[i].len() {
            print!("{} ", matrix[i][j]);
        }
        println!("");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {

        let A: [[f64; 4]; 3] =
            [
                [1., 2., 3., 1.],
                [2., -1., 2., 6.],
                [1., 1., 5., -1.]
            ];

        let result = gauss(&mut A.map(|el| el.to_vec()).to_vec());
        let correct: Vec<f64> = vec!(4.0, -0.0, -1.0);

        assert_eq!(result, correct);
    }
}
