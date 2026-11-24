use std::cell::Cell;

use diffsol::nonlinear_solver::{
    convergence::{Convergence, ConvergenceStatus},
    newton::newton_iteration,
};
use diffsol::{BacktrackingLineSearch, FaerVec, LineSearch, Vector};
use diffsol_nl::{NlError, NonLinearSolverError};

fn vector(value: f64) -> FaerVec<f64> {
    FaerVec::from_vec(vec![value], Default::default())
}

#[test]
fn backtracking_measures_contraction_for_each_affine_solve() {
    let atol = vector(1e-6);
    let mut convergence = Convergence::new(0.0, &atol);
    let mut search = BacktrackingLineSearch::<FaerVec<f64>>::default();
    let mut point = vector(1.0);
    let mut correction = vector(0.0);
    let error_scale = vector(1.0);
    let calls = Cell::new(0);
    for root in [0.0, 2.0] {
        let residual = |x: &FaerVec<f64>, out: &mut FaerVec<f64>| {
            calls.set(calls.get() + 1);
            out[0] = x[0] - root;
        };
        newton_iteration(
            &mut point,
            &mut correction,
            &error_scale,
            residual,
            |_delta| Ok(()),
            &mut convergence,
            &mut search,
        )
        .unwrap();
        assert_eq!(point[0], root);
    }
    // the second solve starts with the rate measured by the first; its first correction is
    // applied and the returned point is still checked with a fresh correction
    assert_eq!(
        calls.get(),
        5,
        "each affine solve must certify its returned point with a fresh correction"
    );
}

#[test]
fn backtracking_keeps_armijo_reduction_for_a_large_nonlinear_step() {
    let atol = vector(1e-10);
    let mut convergence = Convergence::new(0.0, &atol);
    let mut search = BacktrackingLineSearch::<FaerVec<f64>>::default();
    let mut point = vector(0.1);
    let mut correction = vector(0.0);
    let status = search
        .take_optimal_step(
            &mut point,
            &mut correction,
            &vector(1.0),
            &|x, out| out[0] = x[0] * x[0] * x[0] - 1.0,
            &|delta| {
                delta[0] /= 0.03;
                Ok(())
            },
            &mut convergence,
        )
        .unwrap();
    assert!(matches!(status, ConvergenceStatus::Continue));
    assert!(search.n_iters > 0, "the full Newton step must be reduced");
    assert!((point[0].powi(3) - 1.0).abs() < (0.1_f64.powi(3) - 1.0).abs());
}

#[test]
fn backtracking_applies_the_correction_whose_remaining_error_it_certifies() {
    let atol = vector(1e-3);
    let mut convergence = Convergence::with_tolerance(0.0, &atol, 0.2);
    let mut search = BacktrackingLineSearch::<FaerVec<f64>>::default();
    let mut point = vector(1.0);
    let mut correction = vector(0.0);
    let calls = Cell::new(0);
    newton_iteration(
        &mut point,
        &mut correction,
        &vector(1.0),
        |x, out| {
            calls.set(calls.get() + 1);
            out[0] = x[0] + 0.01 * x[0] * x[0];
        },
        |delta| {
            delta[0] /= 1.02;
            Ok(())
        },
        &mut convergence,
        &mut search,
    )
    .unwrap();
    // The exact nearby root is zero. Applying only the first correction gives
    // x=0.0098, even though the next correction's rate estimate is below tol.
    assert!(
        point[0].abs() <= 0.2 * atol[0],
        "accepted point {} is outside the requested nonlinear error bound",
        point[0]
    );
    assert_eq!(
        calls.get(),
        3,
        "the returned point needs one fresh correction check after measured contraction"
    );
}

#[test]
fn retained_rate_does_not_certify_a_different_nonlinear_solve() {
    let atol = vector(1e-6);
    let mut convergence = Convergence::new(0.0, &atol);
    let mut search = BacktrackingLineSearch::<FaerVec<f64>>::default();
    let mut point = vector(1.0);
    let mut correction = vector(0.0);
    let scale = vector(1.0);
    newton_iteration(
        &mut point,
        &mut correction,
        &scale,
        |x, out| out[0] = x[0],
        |_delta| Ok(()),
        &mut convergence,
        &mut search,
    )
    .unwrap();
    assert_eq!(point[0], 0.0);
    point[0] = 1.0;
    newton_iteration(
        &mut point,
        &mut correction,
        &scale,
        |x, out| out[0] = x[0] + 0.01 * x[0] * x[0],
        |_delta| Ok(()),
        &mut convergence,
        &mut search,
    )
    .unwrap();
    assert!(
        point[0].abs() <= 0.2 * atol[0],
        "stale affine rate accepted x={} outside current nonlinear budget",
        point[0]
    );
}

#[test]
fn backtracking_does_not_certify_a_point_after_crossing_into_a_slow_branch() {
    let atol = vector(1e-6);
    let mut convergence = Convergence::new(0.0, &atol);
    let mut search = BacktrackingLineSearch::<FaerVec<f64>>::default();
    let mut point = vector(1.0);
    let mut correction = vector(0.0);
    let boundary = 1e-5 / 0.99;
    let result = newton_iteration(
        &mut point,
        &mut correction,
        &vector(1.0),
        |x, out| {
            out[0] = if x[0] >= boundary {
                x[0]
            } else {
                0.01 * x[0] + 1e-5
            }
        },
        |_delta| Ok(()),
        &mut convergence,
        &mut search,
    );
    // The continuous, strictly increasing function has its unique root at
    // -.001. The cached inverse Jacobian is one, valid for the initial branch.
    // Its .99 contraction after the branch change cannot converge in ten steps.
    if result.is_ok() {
        assert!(
            (point[0] + 0.001).abs() <= 0.2 * atol[0],
            "false convergence at{} after crossing a slope change; exact root-.001",
            point[0]
        );
    } else {
        assert!(matches!(
            result,
            Err(NlError::NonLinearSolverError(
                NonLinearSolverError::NewtonDiverged | NonLinearSolverError::NewtonMaxIterations
            ))
        ));
    }
}
