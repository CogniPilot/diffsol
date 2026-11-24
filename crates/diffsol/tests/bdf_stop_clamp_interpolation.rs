use diffsol::{FaerLU, FaerMat, FaerVec, OdeBuilder, OdeSolverMethod, Vector};

fn values(v: &FaerVec<f64>) -> Vec<f64> {
    (0..v.len()).map(|i| v.get_index(i)).collect()
}

#[test]
fn native_bdf_previous_point_survives_a_close_hard_stop() {
    // CanyonRun.Run, stage-instrumented: an order-5 step was accepted
    // just short of a hard stop, the closing stop clamp shrank the step about
    // 2500x, and the previous point was then extrapolated through ~1e14
    // interpolation weights. Here a stop 1e-4 of a step ahead does the same.
    let problem = OdeBuilder::<FaerMat<f64>>::new()
        .t0(0.0)
        .h0(1e-4)
        .rtol(1e-6)
        .atol([1e-6, 1e-6])
        .rhs_implicit(
            |x, _p, _t, out| {
                out[0] = x[1];
                out[1] = -x[0];
            },
            |_x, _p, _t, v, out| {
                out[0] = v[1];
                out[1] = -v[0];
            },
        )
        .init(
            |_p, _t, out| {
                out[0] = 60.0;
                out[1] = 15.0;
            },
            2,
        )
        .build()
        .unwrap();
    // One clamp per fresh integration: a clamp leaves the next step tiny, and
    // clamping again from there would put the stop inside Diffsol's own
    // roundoff window (StopTimeAtCurrentTime) instead of testing interpolation.
    for extra_steps in 0..5 {
        let mut method = problem.bdf::<FaerLU<f64>>().unwrap();
        let mut reached = None;
        let mut clamp_point = None;
        for _ in 0..400 {
            let previous = (method.state().t, values(method.state().y));
            method.step().unwrap();
            if method.order() < 4 {
                continue;
            }
            let seen = reached.get_or_insert(0);
            if *seen == extra_steps {
                clamp_point = Some(previous);
                break;
            }
            *seen += 1;
        }
        assert_eq!(
            reached,
            Some(extra_steps),
            "the run must reach order four and clamp after {extra_steps} more steps"
        );
        // The step loop stopped right after the chosen accepted step, so the
        // clamp below sees exactly the state it would have seen inside it.
        let (previous_time, previous_states) = clamp_point
            .unwrap_or_else(|| panic!("the run must clamp after {extra_steps} more steps"));
        let (time, step) = (method.state().t, method.state().h);
        method.set_stop_time(time + 1e-4 * step).unwrap();
        let sampled = method.interpolate(previous_time).unwrap();
        for (index, (expected, actual)) in previous_states
            .iter()
            .zip(values(&sampled).iter())
            .enumerate()
        {
            assert!(
                (expected - actual).abs() <= 1e-9 * expected.abs().max(1.0),
                "state {index} at the previous point {previous_time} after a stop \
                 clamp at order {}: accepted={expected}, interpolated={actual}",
                method.order()
            );
        }
    }
}
