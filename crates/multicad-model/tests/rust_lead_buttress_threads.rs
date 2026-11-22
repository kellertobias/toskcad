use multicad_core::{
    commands::{Actor, Channel, Session},
    document::{Document, EntityId},
    model::{
        self, BodyInput, ButtressLoadFlankSide, PrimitiveParameters, PrimitiveRequest,
        ThreadFamily, ThreadHand, ThreadMode, ThreadParameters, ThreadRequest, length,
    },
    project_file::{self, ProjectFile},
};
use multicad_model::{Compiler, GeometryBackend};
use multicad_precision::PortableMath as _;
fn actor() -> Actor {
    Actor::local(Channel::Mcp)
}
fn compile(
    s: &Session,
    backend: GeometryBackend,
) -> Result<multicad_model::Model, multicad_model::ModelError> {
    Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile_with_backend(s.document(), s.revision(), backend, None)
}
fn reopen(s: &Session) -> Session {
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(s.document().clone()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    Session::new(
        project_file::read(
            &mut bytes.as_slice(),
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .unwrap()
        .project
        .document,
    )
    .unwrap()
}

#[test]
fn sourced_lead_and_asymmetric_root_pairs_keep_geometry_phase_history_and_step() {
    let mut cases: Vec<ThreadParameters> = model::lead_screw_presets()
        .into_iter()
        .map(|row| row.parameters())
        .collect();
    let mut multi = model::lead_screw_presets().remove(1).parameters();
    multi.hand = ThreadHand::Left;
    multi.starts = 2;
    cases.push(multi);
    for side in [
        ButtressLoadFlankSide::NegativeAxisSide,
        ButtressLoadFlankSide::PositiveAxisSide,
    ] {
        for hand in [ThreadHand::Right, ThreadHand::Left] {
            let mut p = model::buttress_thread_presets().remove(0).parameters(side);
            p.hand = hand;
            cases.push(p);
        }
    }

    for (index, mut p) in cases.into_iter().enumerate() {
        let mut s =
            Session::new(Document::new("Rust lead and buttress geometry").unwrap()).unwrap();
        if index % 2 != 0 {
            p.frame = model::PrimitiveFrame {
                origin: [17., -8., 13.].map(length),
                u: [0.8, 0.6, 0.],
                v: [-0.36, 0.48, 0.8],
            };
        }
        let initial = p.resolve(s.parameters(), ThreadMode::NewExternal).unwrap();
        p.length = length(3.37 * initial.pitch_mm);
        let major = p.major_diameter.clone();
        let pitch = p.pitch.clone();
        let male = ThreadRequest::new("Sourced external", p.clone());
        s.commit(
            model::prepare_thread(&s, s.revision(), actor(), &male)
                .unwrap()
                .prepared,
        )
        .unwrap();
        let resolved = p.resolve(s.parameters(), ThreadMode::NewExternal).unwrap();
        let blank = PrimitiveRequest::new(
            "Receiving blank",
            PrimitiveParameters::Cylinder {
                radius: length(resolved.major_radius_mm + 3.),
                height: p.length.clone(),
                frame: p.frame.clone(),
            },
        );
        s.commit(
            model::prepare_primitive(&s, s.revision(), actor(), &blank)
                .unwrap()
                .prepared,
        )
        .unwrap();
        let female = ThreadRequest::receiving(
            "Sourced internal",
            ThreadMode::CutInternal,
            BodyInput {
                body: blank.body,
                producer: blank.feature,
            },
            p.clone(),
        )
        .unwrap();
        s.commit(
            model::prepare_thread(&s, s.revision(), actor(), &female)
                .unwrap()
                .prepared,
        )
        .unwrap();
        let rust =
            compile(&s, GeometryBackend::Rust).unwrap_or_else(|e| panic!("case{index}: {e}"));

        for (internal, request) in [(false, &male), (true, &female)] {
            let solid = rust.solid(request.feature).unwrap();
            assert!(solid.uses_rust_geometry());
            let r = request
                .parameters
                .resolve(s.parameters(), request.mode)
                .unwrap();
            let world = |point: [f64; 3]| {
                let vector = r.frame.vector(point);
                std::array::from_fn(|i| r.frame.origin[i] + vector[i])
            };
            let phase = r.buttress.as_ref().map_or(0., |b| b.crest_phase_axial_mm);
            for (sign, angle) in [
                (-1., r.profile.lower_flank_deg),
                (1., r.profile.upper_flank_deg),
            ] {
                for fraction in [0.2, 0.6] {
                    let drop = r.profile.depth_mm * fraction;
                    let z = 2. * r.pitch_mm
                        + phase
                        + sign
                            * (r.profile.crest_width_mm / 2.
                                + drop * angle.to_radians().portable_tan());
                    for (offset, in_profile) in [(-0.001, true), (0.001, false)] {
                        let point = [r.major_radius_mm - drop + offset, 0., z];

                        assert_eq!(
                            solid.rust_thread_contains(world(point)).unwrap(),
                            if internal { !in_profile } else { in_profile },
                            "case{index} internal{internal} flank{sign} fraction{fraction}"
                        );
                    }
                }
            }
            if let Some(circle) = r.buttress.as_ref().and_then(|b| b.root_circle.as_ref()) {
                assert!(!internal);
                for angle in [-0.5_f64, 0.5] {
                    let z = 2. * r.pitch_mm
                        + circle.center_axial_mm
                        + circle.radius_mm * angle.portable_sin();
                    let boundary =
                        circle.center_radius_mm - circle.radius_mm * angle.portable_cos();

                    assert!(
                        solid
                            .rust_thread_contains(world([boundary - 0.001, 0., z]))
                            .unwrap()
                    );
                    assert!(
                        !solid
                            .rust_thread_contains(world([boundary + 0.001, 0., z]))
                            .unwrap()
                    );
                }
            }
            let folder = tempfile::tempdir().unwrap();
            let path = folder.path().join("thread.step");
            solid
                .write_step(&path)
                .unwrap_or_else(|e| panic!("STEP case{index} internal{internal}: {e}"));
            let verified = solid
                .verify_helical_step(&std::fs::read_to_string(&path).unwrap())
                .unwrap()
                .unwrap();
            assert!((verified - solid.volume_mm3().unwrap()).abs() < 1e-3);
            let imported = multicad_kernel::Solid::import_rust_step(
                &std::fs::read_to_string(&path).unwrap(),
                EntityId::new(),
                Default::default(),
            )
            .unwrap();
            assert!((imported.volume_mm3().unwrap() - solid.volume_mm3().unwrap()).abs() < 1e-3);
        }

        let before = s.document().clone();
        let updated = ThreadParameters {
            length: length(3.75 * resolved.pitch_mm),
            hand: if p.hand == ThreadHand::Right {
                ThreadHand::Left
            } else {
                ThreadHand::Right
            },
            ..p.clone()
        };
        s.commit(
            model::prepare_thread_edit(&s, s.revision(), actor(), male.feature, &updated)
                .unwrap()
                .prepared,
        )
        .unwrap();
        let after = s.document().clone();
        let opened = reopen(&s);
        let saved = model::read_thread(opened.document(), male.feature).unwrap();
        assert_eq!(saved.parameters, updated);
        assert_eq!(saved.parameters.major_diameter, major);
        assert_eq!(saved.parameters.pitch, pitch);
        let rebuilt = compile(&opened, GeometryBackend::Rust).unwrap();
        assert!(
            (rebuilt.solid(male.feature).unwrap().volume_mm3().unwrap()
                - rust.solid(male.feature).unwrap().volume_mm3().unwrap() * 3.75 / 3.37)
                .abs()
                < 1e-6
        );
        s.undo(s.revision(), actor()).unwrap();
        assert_eq!(s.document(), &before);
        s.redo(s.revision(), actor()).unwrap();
        assert_eq!(s.document(), &after);
        let invalid = ThreadParameters {
            minimum_wall: length(10.),
            ..p
        };
        let candidate =
            model::prepare_thread_edit(&s, s.revision(), actor(), female.feature, &invalid)
                .unwrap();
        assert!(
            Compiler {
                omit_display: true,
                ..Default::default()
            }
            .compile_with_backend(
                candidate.prepared.document(),
                s.revision() + 1,
                GeometryBackend::Rust,
                None
            )
            .is_err()
        );
        assert_eq!(s.document(), &after);
    }
}

#[test]
fn sourced_rounded_and_tapered_families_keep_their_distinct_admitted_laws() {
    let mut cases: Vec<_> = model::whitworth_thread_presets()
        .into_iter()
        .map(|row| ThreadParameters {
            family: row.family,
            fit: if row.family == ThreadFamily::Bspp {
                model::ThreadFit::BsppClassBExternal
            } else {
                model::ThreadFit::WhitworthMediumExternal
            },
            major_diameter: length(row.diameter_mm),
            pitch: length(row.pitch_mm),
            length: length(3. * row.pitch_mm),
            ..Default::default()
        })
        .collect();
    cases.push(
        model::PipeThreadIntent {
            form: model::PipeThreadForm::R,
            size: model::PipeThreadSize::OneEighth,
            gauge: model::PipeGaugeConvention::Iso7Nominal1994,
        }
        .parameters(),
    );
    cases.push(
        model::NptThreadIntent {
            form: model::NptThreadForm::External,
            size: model::NptThreadSize::OneEighth,
            gauge: model::NptGaugeConvention::Nbs1959PlainFace,
            truncation: model::NptTruncationTarget::MidpointConeParallel,
        }
        .parameters(),
    );
    for parameters in cases {
        let s = Session::new(Document::new("Distinct rounded or tapered family").unwrap()).unwrap();
        let family = parameters.family;
        let resolved = parameters
            .resolve(s.parameters(), ThreadMode::NewExternal)
            .unwrap();
        let request = ThreadRequest::new("Sourced profile", parameters);
        let candidate = model::prepare_thread(&s, 0, actor(), &request).unwrap();
        let compiled = Compiler {
            omit_display: true,
            ..Default::default()
        }
        .compile_with_backend(
            candidate.prepared.document(),
            1,
            GeometryBackend::Rust,
            None,
        );
        if family == ThreadFamily::Npt {
            let compiled = compiled.unwrap();
            let solid = compiled.solid(request.feature).unwrap();
            let clipped = solid
                .rust_clipped_npt_source()
                .expect("Actual clipped material, not the untrimmed pipe provider");
            let spec = clipped.source().spec();
            assert_eq!(spec.major_radius_mm, resolved.major_radius_mm);
            assert_eq!(spec.pitch_mm, resolved.pitch_mm);
            assert_eq!(spec.length_mm, resolved.length_mm);
            let limit = resolved
                .npt
                .as_ref()
                .unwrap()
                .stock_crest_radius_limit_mm
                .unwrap();
            assert_eq!(clipped.stock_radius_limit_mm(), limit);
            // Four thread bands, the actual stock cylinder and two end caps.
            // The independent STEP rectangular partition has separate owners.
            assert_eq!(solid.faces().unwrap().len(), 7);
            assert!(solid.rust_clipped_npt_face_measures().unwrap().iter().any(|face| face.face() == multicad_kernel::rust_kernel::axial_pipe_thread::authored::StockClipFace::StockCylinder));
            assert!(
                solid
                    .rust_clipped_npt_contains([0., 0., resolved.length_mm / 2.])
                    .unwrap()
            );
            assert!(
                !solid
                    .rust_clipped_npt_contains([limit + 0.01, 0., resolved.length_mm / 2.])
                    .unwrap()
            );
            let bounds = solid.bounds().unwrap();
            assert!((bounds.max[0] - limit).abs() < 1e-10);
            assert!((bounds.max[2] - resolved.length_mm).abs() < 1e-10);
            let [lo, hi] = solid.rust_clipped_npt_volume_bounds_mm3().unwrap();
            let volume = solid.volume_mm3().unwrap();
            assert!(lo > 0. && lo <= volume && volume <= hi);
        } else {
            let compiled = compiled.unwrap();
            let solid = compiled.solid(request.feature).unwrap();
            assert!(solid.uses_rust_geometry());
            if family == ThreadFamily::Bspt {
                let (source, _) = solid
                    .rust_authored_axial_pipe_thread_source()
                    .expect("Actual authored pipe law");
                let spec = source.transport.spec();
                assert_eq!(spec.major_radius_mm, resolved.major_radius_mm);
                assert_eq!(spec.pitch_mm, resolved.pitch_mm);
                assert_eq!(spec.length_mm, resolved.length_mm);
                assert_eq!(spec.radial_taper, resolved.profile.radial_taper);
                assert!(matches!(spec.section, multicad_kernel::rust_kernel::axial_pipe_thread::authored::AuthoredSection::BsptRounded { .. }));
            } else {
                let data = solid
                    .rust_rounded_thread_data()
                    .expect("Actual two-circle rounded law");
                let spec = data.thread.spec();
                assert_eq!(spec.major_radius_mm, resolved.major_radius_mm);
                assert_eq!(spec.pitch_mm, resolved.pitch_mm);
                assert_eq!(spec.length_mm, resolved.length_mm);
                assert_eq!(spec.depth_mm, resolved.profile.depth_mm);
                assert_eq!(
                    Some(spec.crest_radius_mm),
                    resolved.profile.rounding_radius_mm
                );
                assert_eq!(spec.root_radius_mm, spec.crest_radius_mm);
                assert_eq!(spec.flank_angle_deg, resolved.profile.lower_flank_deg);
            }
            assert!(solid.volume_mm3().unwrap() > 0.);
        }
        assert_eq!(s.revision(), 0);
    }
}
