use crate::field::Field;

pub fn empty(_: &mut Field, _: &Field, _: &Field, _: &Field, _: f64, _: f64, _: &Fn(f64, f64, &Field) -> f64, _: &Fn(f64, f64, &Fn(f64, f64) -> f64, f64) -> f64) {

}

pub fn semi_lagrangian(field_dst: &mut Field, field_src: &Field, x_velocity: &Field, y_velocity: &Field, dt: f64, dx: f64, interpolation: &Fn(f64, f64, &Field) -> f64, integration: &Fn(f64, f64, &Fn(f64, f64) -> f64, f64) -> f64) {
    for row in 0..field_dst.rows {
        for column in 0..field_dst.columns {

            // get position on grid
            let x = column as f64 + field_dst.offset_x;
            let y = row as f64 + field_dst.offset_y;

//            let f1 = |_: f64, o: f64| -interpolation(o, x, &x_velocity) / dx;
//            let f2 = |_: f64, o: f64| -interpolation(y, o, &y_velocity) / dx;

             // set up function which interpolates
            let f1 = |_t: f64, y0: f64| -1.0 * interpolation(y0, y, &x_velocity) / dx;
            let f2 = |_t: f64, y0: f64| -1.0 * interpolation(x, y0, &y_velocity) / dx;

            let old_x: f64 = integration(0.0, x, &f1, dt);
            let old_y: f64 = integration(0.0, y, &f2, dt);

//            println!("x: {}, y: {}, oldx: {}, oldy: {}", x, y, old_x, old_y);


            *field_dst.at_mut(row, column) = interpolation(old_x, old_y, field_src);

//            runge_kutta3(x_velocity, y_velocity, &x, &y, dt, dx);
//
//            *field.at_mut(row, column) = lerp_grid()

        }
    }
}

//fn runge_kutta3(u: &Field, v: &Field, x: &mut f64, y: &mut f64, time_step: f64, dx: f64) {
//    let x_velocity1 = lerp_grid(&u, *x, *y) / dx;
//    let y_velocity1 = lerp_grid(&v, *x, *y) / dx;
//
//    let mid_x = *x - 0.5 * time_step * x_velocity1;
//    let mid_y = *y - 0.5 * time_step * y_velocity1;
//
//    let x_velocity2 = lerp_grid(&u, mid_x, mid_y) / dx;
//    let y_velocity2 = lerp_grid(&v, mid_x, mid_y) / dx;
//
//    let last_x = *x - 0.75 * time_step * x_velocity2;
//    let last_y = *y - 0.75 * time_step * y_velocity2;
//
//    let x_velocity3 = lerp_grid(&u, last_x, last_y) / dx;
//    let y_velocity3 = lerp_grid(&v, last_x, last_y) / dx;
//
//    *x -= time_step * ((2.0/9.0) * x_velocity1 + (3.0/9.0) * x_velocity2 + (4.0/9.0) * x_velocity3);
//    *y -= time_step * ((2.0/9.0) * y_velocity1 + (3.0/9.0) * y_velocity2 + (4.0/9.0) * y_velocity3);
//}
//
//pub fn lerp_grid(grid: &Field, x: f64, y: f64) -> f64 {
//    let mut temp_x = min(max(x - grid.offset_x, 0.0), grid.columns as f64 - 1.001);
//    let mut temp_y = min(max(y - grid.offset_y, 0.0), grid.rows as f64 - 1.001);
//
//    let ix = temp_x as usize;
//    let iy = temp_y as usize;
//
//    temp_x -= ix as f64;
//    temp_y -= iy as f64;
//
//    let x00: f64 = grid.at(iy + 0, ix + 0);
//    let x10: f64 = grid.at(iy + 0, ix + 1);
//    let x01: f64 = grid.at(iy + 1, ix + 0);
//    let x11: f64 = grid.at(iy + 1, ix + 1);
//
//    lerp(lerp(x00, x10, temp_x), lerp(x01, x11, temp_x), temp_y)
//}
//
//pub fn lerp(a: f64, b: f64, c: f64) -> f64 {
//    a * (1.0 - c) + b * c
//}
//
//pub fn min<T: PartialOrd>(a: T, b: T) -> T {
//    if a < b {
//        return a
//    }
//    return b
//}
//
//pub fn max<T: PartialOrd>(a: T, b: T) -> T {
//    if a > b {
//        return a
//    }
//    return b
//}