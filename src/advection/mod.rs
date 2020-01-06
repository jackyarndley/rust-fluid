use crate::interpolation::Interpolation;
use crate::integration::Integration;
use crate::util::fluid_quantity::{FluidQuantity};

pub enum Advection {
    SemiLagrangian
}

impl Advection {
    pub fn advect(&self, u_velocity: &mut FluidQuantity, v_velocity: &mut FluidQuantity, density: &mut FluidQuantity, timestep: f64, interpolation: &Interpolation, integration: &Integration) {
        match self {
            Advection::SemiLagrangian => {
                for row in 0..u_velocity.rows {
                    for column in 0..u_velocity.columns {
                        let x = column as f64 + u_velocity.x_offset;
                        let y = row as f64 + u_velocity.y_offset;

                        let f1 = |_t: f64, y0: f64| -1.0 * interpolation.run(y0, y, &u_velocity) / timestep;
                        let f2 = |_t: f64, y0: f64| -1.0 * interpolation.run(x, y0, &v_velocity) / timestep;

                        let old_x: f64 = integration.run(0.0, x, &f1, timestep);
                        let old_y: f64 = integration.run(0.0, y, &f2, timestep);

                        *u_velocity.dst_at_mut(row, column) = interpolation.run(old_x, old_y, u_velocity);
                    }
                }

                for row in 0..v_velocity.rows {
                    for column in 0..v_velocity.columns {
                        let x = column as f64 + v_velocity.x_offset;
                        let y = row as f64 + v_velocity.y_offset;

                        let f1 = |_t: f64, y0: f64| -1.0 * interpolation.run(y0, y, &u_velocity) / timestep;
                        let f2 = |_t: f64, y0: f64| -1.0 * interpolation.run(x, y0, &v_velocity) / timestep;

                        let old_x: f64 = integration.run(0.0, x, &f1, timestep);
                        let old_y: f64 = integration.run(0.0, y, &f2, timestep);

                        *v_velocity.dst_at_mut(row, column) = interpolation.run(old_x, old_y, v_velocity);
                    }
                }

                for row in 0..density.rows {
                    for column in 0..density.columns {
                        let x = column as f64 + density.x_offset;
                        let y = row as f64 + density.y_offset;

                        let f1 = |_t: f64, y0: f64| -1.0 * interpolation.run(y0, y, &u_velocity) / timestep;
                        let f2 = |_t: f64, y0: f64| -1.0 * interpolation.run(x, y0, &v_velocity) / timestep;

                        let old_x: f64 = integration.run(0.0, x, &f1, timestep);
                        let old_y: f64 = integration.run(0.0, y, &f2, timestep);

                        *density.dst_at_mut(row, column) = interpolation.run(old_x, old_y, density);
                    }
                }

                u_velocity.swap_buffers();
                v_velocity.swap_buffers();
                density.swap_buffers();
            }
        }
    }
}