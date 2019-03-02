use rust_fluid::fluid_solver::FluidSolver;
use rust_fluid::{interpolation, integration, linear_solvers, advection};

fn main() {
    let mut buffer = vec![0u8; 128 * 128 * 4];

    let mut solver = FluidSolver::new(128, 128, 0.005, 1.0 / 128.0, 1.0)
        .interpolation(interpolation::bicubic_interpolate)
        .integration(integration::bogacki_shampine)
        .linear_solver(linear_solvers::conjugate_gradient)
        .advection(advection::semi_lagrangian);

    solver.u_velocity.add_inflow(20, 20, 40, 40, 2.0);
    solver.v_velocity.add_inflow(20, 20, 40, 40, 2.0);
    solver.density.add_inflow(20, 20, 40, 40, 2.0);

    for iteration in 0..1000 {
        for i in 0..4 {
            print!("Global iteration {}. ", 4 * iteration + i);
            solver.solve();


            solver.u_velocity.add_inflow(20, 20, 40, 40, 2.0);
            solver.v_velocity.add_inflow(20, 20, 40, 40, 2.0);
            solver.density.add_inflow(20, 20, 40, 40, 2.0);
        }

        solver.to_image(&mut buffer);
        lodepng::encode32_file(format!("output/{}.png", iteration), &buffer, 128, 128).unwrap();
    }
}
