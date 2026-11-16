use super::*;
use crate::{
    axial_pipe_thread::authored::StockClipPlacement,
    orthogonal::{self, Box3, Operation},
    orthogonal_bore_seats::{Options as SourceOptions, Seat, Seats},
    orthogonal_boundary, orthogonal_pockets,
    step::orthogonal_seats,
};
use std::{cell::Cell, rc::Rc};
fn text() -> String {
    text_placed(
        StockClipPlacement::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0.; 3]).unwrap(),
    )
}
pub(super) fn text_placed(placement: StockClipPlacement) -> String {
    let source = source();
    orthogonal_seats::candidate(
        &source,
        placement,
        orthogonal_seats::Options::default(),
        &|| false,
    )
    .unwrap()
    .as_str()
    .to_owned()
}
pub(super) fn source() -> Seats {
    let stock = Rc::new(
        orthogonal::boolean(
            &[
                Box3 {
                    min: [0.; 3],
                    max: [40., 40., 5.],
                },
                Box3 {
                    min: [0.; 3],
                    max: [40., 5., 40.],
                },
            ],
            &[],
            Operation::Union,
            1e-7,
        )
        .unwrap()
        .unwrap(),
    );
    let boundary = orthogonal_boundary::build(stock.quads()).unwrap();
    let bores = Rc::from(
        [
            [[20., 20., 5.1], [20., 20., -0.1]],
            [[20., 5.1, 20.], [20., -0.1, 20.]],
        ]
        .map(|[a, b]| {
            orthogonal_pockets::identify_through_bore(&stock, &boundary, a, b, 4.9, 1e-7, 1e-5)
                .unwrap()
        }),
    );
    let seats = [0, 1].map(|bore| Seat {
        bore,
        end: false,
        radius_mm: 8.02,
        depth_mm: 3.12,
    });
    Seats::new(
        stock,
        bores,
        &seats,
        SourceOptions {
            clearance_mm: 1e-7,
            minimum_feature_mm: 1e-5,
            maximum_work: 16_000_000,
            maximum_bytes: 64 * 1024 * 1024,
        },
        &|| false,
    )
    .unwrap()
}
fn parse_document(text: &str) -> CompoundDocument<'_> {
    CompoundDocument::parse(text, owned_insert::DocumentOptions::default(), None).unwrap()
}
#[test]
fn literal_perpendicular_seats_keep_original_owners_and_decimal_spellings() {
    let text = text();
    let document = parse_document(&text);
    let item = document.represented_members()[0].item_entity();
    let boundary = decode(&document, item, Options::default(), &|| false).unwrap();
    assert!(std::ptr::eq(boundary.document(), &document));
    assert_eq!(boundary.item_entity(), item);
    assert_eq!(boundary.solid_entity(), item);
    assert_eq!(boundary.planar_face_count(), 8);
    assert_eq!(boundary.conical_face_count(), 2);
    assert_eq!(boundary.cylindrical_face_count(), 2);
    assert_eq!(boundary.face_entities().len(), 12);
    let chains = boundary
        .through_chains(Options::default(), &|| false)
        .unwrap();
    assert!(std::ptr::eq(chains.boundary(), &boundary));
    assert_eq!(chains.bores().len(), 2);
    let mut owned = BTreeSet::new();
    for bore in chains.bores() {
        assert_eq!(bore.wall_face_entities().len(), 2);
        assert_eq!(bore.rim_edge_entities().len(), 3);
        assert_eq!(bore.rim_circles().len(), 3);
        assert!(bore.cylinder_index() < 2);
        assert_ne!(bore.cap_face_entities()[0], bore.cap_face_entities()[1]);
        for face in bore.wall_face_entities() {
            assert!(owned.insert(*face));
        }
    }
    assert_eq!(owned.len(), 4);
    assert!(boundary.edge_count() > 20 && boundary.vertex_count() > 10);
    assert!(
        boundary
            .numeric_leaves()
            .iter()
            .all(|n| boundary.dependencies().binary_search(&n.entity()).is_ok())
    );
    let leaf = boundary
        .numeric_leaves()
        .iter()
        .find(|n| n.spelling().contains("4.9000000000000004"))
        .unwrap();
    let record = &document.records()[&leaf.entity()][0];
    let spellings = part21::numeric_spellings(record).unwrap();
    assert_eq!(spellings[leaf.ordinal()], leaf.spelling());
    let original = leaf.spelling();
    // Append a decimal mantissa digit, preserving the parsed f64 cache while
    // changing the actual literal. This stage retains both, never equates them.
    let (mantissa, exponent) = original.split_once('E').unwrap();
    let altered = format!("{mantissa}0000000000000001E{exponent}");
    assert_eq!(
        original.parse::<f64>().unwrap(),
        altered.parse::<f64>().unwrap()
    );
    let changed = text.replacen(original, &altered, 1);
    let changed_document = parse_document(&changed);
    let changed_boundary = decode(&changed_document, item, Options::default(), &|| false).unwrap();
    assert!(
        changed_boundary
            .numeric_leaves()
            .iter()
            .any(|n| n.spelling() == altered)
    );
}
#[test]
fn malformed_finite_orientation_and_actual_angular_units_are_rejected() {
    let original = text();
    let face = original
        .lines()
        .find(|l| l.contains("ADVANCED_FACE(") && l.ends_with(".T.);"))
        .unwrap();
    let flipped = face.replacen(".T.)", ".F.)", 1);
    assert_ne!(flipped, face);
    let changed = original.replacen(face, &flipped, 1);
    let document = parse_document(&changed);
    let item = document.represented_members()[0].item_entity();
    assert!(matches!(
        decode(&document, item, Options::default(), &|| false),
        Err(Error::Invalid)
    ));
    let degrees = original.replace(".RADIAN.", ".DEGREE.");
    assert_ne!(degrees, original);
    let document = parse_document(&degrees);
    let item = document.represented_members()[0].item_entity();
    assert!(matches!(
        decode(&document, item, Options::default(), &|| false),
        Err(Error::Unsupported)
    ));
}
#[test]
fn selected_query_has_separate_exact_work_frontier_and_atomic_cancellation() {
    let original = text();
    let document = parse_document(&original);
    let item = document.represented_members()[0].item_entity();
    let earned = decode(&document, item, Options::default(), &|| false).unwrap();
    let exact = Options {
        maximum_work: earned.work_used(),
        maximum_bytes: earned.query_budget_bytes(),
    };
    assert!(decode(&document, item, exact, &|| false).is_ok());
    assert!(matches!(
        decode(
            &document,
            item,
            Options {
                maximum_work: exact.maximum_work - 1,
                ..exact
            },
            &|| false
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        decode(
            &document,
            item,
            Options {
                maximum_bytes: exact.maximum_bytes - 1,
                ..exact
            },
            &|| false
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        decode(
            &document,
            item,
            Options {
                maximum_work: 0,
                ..exact
            },
            &|| false
        ),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        decode(&document, u32::MAX, exact, &|| false),
        Err(Error::Invalid)
    ));
    let chains = earned
        .through_chains(Options::default(), &|| false)
        .unwrap();
    let chain_exact = Options {
        maximum_work: chains.work_used(),
        maximum_bytes: chains.query_budget_bytes(),
    };
    assert_eq!(
        chains.aggregate_known_work(),
        earned.work_used() + chains.work_used()
    );
    assert!(earned.through_chains(chain_exact, &|| false).is_ok());
    assert!(matches!(
        earned.through_chains(
            Options {
                maximum_work: chain_exact.maximum_work - 1,
                ..chain_exact
            },
            &|| false
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        earned.through_chains(
            Options {
                maximum_bytes: chain_exact.maximum_bytes - 1,
                ..chain_exact
            },
            &|| false
        ),
        Err(Error::Limit)
    ));
    let chain_calls = Cell::new(0);
    let chain_cancel = || {
        chain_calls.set(chain_calls.get() + 1);
        chain_calls.get() > 20
    };
    assert!(matches!(
        earned.through_chains(chain_exact, &chain_cancel),
        Err(Error::Cancelled)
    ));
    assert_eq!(chain_calls.get(), 21);
    let calls = Cell::new(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() > 30
    };
    assert!(matches!(
        decode(&document, item, Options::default(), &cancel),
        Err(Error::Cancelled)
    ));
    assert_eq!(calls.get(), 31);
    assert_eq!(
        decode(&document, item, exact, &|| false)
            .unwrap()
            .face_entities(),
        earned.face_entities()
    );
}

#[test]
fn independently_parsed_trim_owner_aliases_and_unsupported_supports_refuse() {
    let original = text();
    let face = original
        .lines()
        .find(|l| l.contains("ADVANCED_FACE("))
        .unwrap();
    let start = face.find(",(").unwrap() + 2;
    let end = start + face[start..].find(')').unwrap();
    let first = face[start..end].split(',').next().unwrap();
    let duplicate = format!(
        "{}{},{}{}",
        &face[..start],
        first,
        &face[start..end],
        &face[end..]
    );
    let duplicated_bound = original.replacen(face, &duplicate, 1);

    let bounds = original
        .lines()
        .filter(|l| l.contains("=FACE_OUTER_BOUND("))
        .take(2)
        .collect::<Vec<_>>();
    let wire = |l: &str| {
        l.split_once(",#")
            .unwrap()
            .1
            .split(',')
            .next()
            .unwrap()
            .to_owned()
    };
    let shared_wire = bounds[1].replacen(
        &format!(",#{},", wire(bounds[1])),
        &format!(",#{},", wire(bounds[0])),
        1,
    );
    let aliased_wire = original.replacen(bounds[1], &shared_wire, 1);
    let unknown_support = original.replacen("CONICAL_SURFACE(", "TOROIDAL_SURFACE(", 1);
    for (label, text, expected) in [
        (
            "duplicate face-bound owner",
            duplicated_bound,
            Error::Invalid,
        ),
        ("shared wire owner", aliased_wire, Error::Invalid),
        ("unsupported support", unknown_support, Error::Unsupported),
    ] {
        assert_ne!(text, original, "{label}");
        let document = parse_document(&text);
        let item = document.represented_members()[0].item_entity();
        assert_eq!(
            decode(&document, item, Options::default(), &|| false).err(),
            Some(expected),
            "{label}"
        );
        // Rejection has no effect on the independently borrowed original file.
        let original_document = parse_document(&original);
        assert!(
            decode(
                &original_document,
                original_document.represented_members()[0].item_entity(),
                Options::default(),
                &|| false
            )
            .is_ok()
        );
    }
}
