use std::mem::swap;

pub struct FaceCenteredVectorGrid2 {
    pub row_data:    Vec<f64>,
    pub column_data: Vec<f64>,
    pub rows:        usize,
    pub columns:     usize,
    pub x_center:    f64,
    pub y_center:    f64
}

impl FaceCenteredVectorGrid2 {
    pub fn new(rows: usize, columns: usize) -> FaceCenteredVectorGrid2 {
        FaceCenteredVectorGrid2 {
            row_data: vec![0.0; rows * (columns + 1)],
            column_data: vec![0.0; (rows + 1) * columns],
            rows,
            columns,
            x_center: 0.5,
            y_center: 0.5
        }
    }

}