use pyo3::prelude::*;

fn kernel_lcs(x: &str, y: &str) -> f64 {
    let x_char: Vec<char> = x.chars().collect();
    let y_char: Vec<char> = y.chars().collect();
    let mut dp = vec![vec![0; y_char.len()+1]; x_char.len()+1];
    for i in 0..x_char.len() {
        for j in 0..y_char.len() {
            if x_char[i] == y_char[j] {
                dp[i+1][j+1] = dp[i][j] + 1;
            } else {
                dp[i+1][j+1] = std::cmp::max(dp[i][j+1], dp[i+1][j]);
            }
        }
    }
    2.0 * dp[x_char.len()][y_char.len()] as f64 / (x_char.len() + y_char.len()) as f64
}

#[pyfunction]
fn gram_matrix_lcs(a: Vec<String>, b: Vec<String>) -> PyResult<Vec<Vec<f64>>> {
    let mut gram_matrix = vec![vec![0.; b.len()]; a.len()];
    for i in 0..a.len() {
        for j in 0..b.len() {
            gram_matrix[i][j] = kernel_lcs(&a[i], &b[j]);
        }
    }
    Ok(gram_matrix)
}

fn kernel_pspectrum(p: usize, x: &str, y: &str) -> f64 {
    let x_char: Vec<char> = x.chars().collect();
    let y_char: Vec<char> = y.chars().collect();
    if x_char.len() < p || y_char.len() < p {
        return -1.0
    }

    let mut memo = std::collections::HashMap::new();
    let mut ret = 0;
    for i in 0..(x_char.len()-p+1) {
        let m: String = x_char[i..i+p].iter().collect();
        let entry = memo.entry(m).or_insert(0);
        *entry += 1;
    }

    for i in 0..(y_char.len()-p+1) {
        let m: String = y_char[i..i+p].iter().collect();
        ret += memo.get(&m).unwrap_or(&0);
    }

    ret as f64 / (((x_char.len() - p + 1) * (y_char.len() - p + 1)) as f64).sqrt()
}

#[pyfunction]
fn gram_matrix_spectrum(p: usize, a: Vec<String>, b: Vec<String>) -> PyResult<Vec<Vec<f64>>> {
    let mut gram_matrix = vec![vec![0.; b.len()]; a.len()];
    for i in 0..a.len() {
        for j in 0..b.len() {
            gram_matrix[i][j] = kernel_pspectrum(p, &a[i], &b[j]);
        }
    }
    Ok(gram_matrix)
}

#[pymodule]
fn str_kernel(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(gram_matrix_lcs, m)?)?;
    m.add_function(wrap_pyfunction!(gram_matrix_spectrum, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::{kernel_lcs, kernel_pspectrum};
    #[test]
    fn test_LCS_kernel() {
        assert_eq!(kernel_lcs("hello", "world"), 0.2);
        assert_eq!(kernel_lcs("hello", "helol"), 0.8);
    }

    #[test]
    fn test_spectrum_kernel() {
        assert_eq!(kernel_pspectrum(2, "hello", "worlld"), 1. / (4.0 * 5.0 as f64).sqrt());
        assert_eq!(kernel_pspectrum(2, "hello", "worlld, hel"), 3. / (4.0 * 10.0 as f64).sqrt());
        assert_eq!(kernel_pspectrum(2, "hello", "hello"), 4.0 / (4.0 * 4.0 as f64).sqrt());
        assert_eq!(kernel_pspectrum(3, "hello, world", "hello, world"), 10.0 / (10.0 * 10.0 as f64).sqrt());
    }
}
