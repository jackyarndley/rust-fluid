use crate::interpolation::Interpolation;
use crate::integration::Integration;
use crate::util::field::Field;

pub enum Advection {
    SemiLagrangian
}

impl Advection {
    pub fn advect(&self, field_dst: &mut Field, field_src: &Field, u_velocity: &Field, v_velocity: &Field, timestep: f64, interpolation: &Interpolation, integration: &Integration) {
        match self {
            Advection::SemiLagrangian => {
                for row in 0..field_dst.rows {
                    for column in 0..field_dst.columns {
                        let x = column as f64 + field_dst.x_offset;
                        let y = row as f64 + field_dst.y_offset;

                        let f1 = |_t: f64, y0: f64| -1.0 * interpolation.run(y0, y, &u_velocity) / timestep;
                        let f2 = |_t: f64, y0: f64| -1.0 * interpolation.run(x, y0, &v_velocity) / timestep;

                        let old_x: f64 = integration.run(0.0, x, &f1, timestep);
                        let old_y: f64 = integration.run(0.0, y, &f2, timestep);

                        *field_dst.at_mut(row, column) = interpolation.run(old_x, old_y, field_src);
                    }
                }
            }
        }
    }
}