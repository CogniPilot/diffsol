use std::cell::Cell;

use diffsol::nonlinear_solver::{
    convergence::{Convergence, ConvergenceStatus},
    newton::newton_iteration,
};
use diffsol::{BacktrackingLineSearch, DiffsolError, FaerVec, LineSearch, Vector};

fn vector(value: f64) -> FaerVec<f64> {
    FaerVec::from_vec(vec![value], Default::default())
}

#[test]
fn backtracking_retains_measured_contraction_between_affine_solves() {
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
    assert_eq!(
        calls.get(),
        3,
        "the second exact affine solve should use its measured contraction history"
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
fn backtracking_rejects_an_iteration_rate_that_cannot_meet_the_budget() {
    let atol = vector(1e-6);
    let mut convergence = Convergence::new(0.0, &atol);
    let mut search = BacktrackingLineSearch::<FaerVec<f64>>::default();
    let mut point = vector(1.0);
    let mut correction = vector(0.0);
    let calls = Cell::new(0);
    let error = newton_iteration(
        &mut point,
        &mut correction,
        &vector(1.0),
        |x, out| {
            calls.set(calls.get() + 1);
            out[0] = x[0];
        },
        |delta| {
            delta[0] /= 20.0;
            Ok(())
        },
        &mut convergence,
        &mut search,
    )
    .unwrap_err();
    assert!(
        matches!(
            error,
            DiffsolError::NonLinearSolverError(
                diffsol::error::NonLinearSolverError::NewtonDiverged
            )
        ),
        "a 0.95 contraction cannot reduce this correction to tolerance in ten iterations"
    );
    assert_eq!(calls.get(), 2);
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
        2,
        "measured fast contraction needs no additional residual calls"
    );
}

/// A backtracking solve whose corrections fit the tolerance but do not
/// contract must not be accepted on a rate retained from an earlier solve:
/// the remaining Newton error is then unbounded by this solve's evidence.
/// (A linear solver that overshoots by 2x makes the iterates alternate about
/// the root with a constant correction, contraction rate 1.)
#[test]
fn backtracking_refuses_corrections_within_tol_that_do_not_contract() {
    let atol = vector(1.0);
    let mut convergence = Convergence::with_tolerance(0.0, &atol, 0.2);
    let mut search = BacktrackingLineSearch::<FaerVec<f64>>::default();
    let mut correction = vector(0.0);
    // An exact solve leaves a small retained rate.
    let mut point = vector(1.0);
    newton_iteration(
        &mut point,
        &mut correction,
        &vector(1.0),
        |x, out| out[0] = x[0],
        |_delta| Ok(()),
        &mut convergence,
        &mut search,
    )
    .unwrap();
    // A non-contracting solve that starts within tolerance.
    let mut point = vector(0.05);
    let outcome = newton_iteration(
        &mut point,
        &mut correction,
        &vector(1.0),
        |x, out| out[0] = x[0],
        |delta| {
            delta[0] *= 2.0;
            Ok(())
        },
        &mut convergence,
        &mut search,
    );
    assert!(
        outcome.is_err(),
        "accepted x = {} although its corrections never contracted (root 0)",
        point[0]
    );
}

/// A first correction of exactly zero means the point already solves the
/// system; it is accepted without a contraction rate (which would be 0/0).
#[test]
fn backtracking_accepts_an_exact_initial_point() {
    let atol = vector(1e-6);
    let mut convergence = Convergence::new(0.0, &atol);
    let mut search = BacktrackingLineSearch::<FaerVec<f64>>::default();
    let mut point = vector(2.0);
    let mut correction = vector(0.0);
    newton_iteration(
        &mut point,
        &mut correction,
        &vector(1.0),
        |x, out| out[0] = x[0] - 2.0,
        |_delta| Ok(()),
        &mut convergence,
        &mut search,
    )
    .unwrap();
    assert_eq!(point[0], 2.0);
}
