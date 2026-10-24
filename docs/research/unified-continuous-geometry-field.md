# Unified Continuous Geometry Field: Research Design

Status: proposal for the first controlled experiment  
Issue: https://github.com/cxs503/TensorCAD/issues/6

## 1. Research decision

TensorCAD should investigate **one representation for every shape**: a continuous, state-conditioned implicit field.

The production representation is a function

\[
F_\theta(\mathbf{x},\mathbf{c}) \rightarrow s
\]

where:

- \(\mathbf{x}=(x,y,z)\) is a point in a canonical shape coordinate system;
- \(\mathbf{c}\) is the editable geometry state (named design parameters and, when needed, compact shape/deformation controls);
- \(F_\theta\) is one shared field model family;
- the zero level set \(F_\theta(\mathbf{x},\mathbf{c})=0\) defines the boundary; the sign defines inside/outside.

The target is a **continuous geometric field**, not a mesh and not a collection of unrelated shape-specific backends. Meshes are sampled outputs for visualization, CAD exchange, and solvers—not the canonical shape representation.

This is a research hypothesis, not a claim that a neural field is already the best answer. The experiment must test whether a single learned field can be accurate, stable, and genuinely editable enough for engineering use.

## 2. One representation, not a dispatcher

Do not implement a rule such as “use an analytic formula for simple shapes, a spline for another shape, and a neural network when formulas fail.” That would preserve several competing representations and merely hide a dispatcher behind an API.

Instead:

1. Every evaluated shape is represented by the same field contract and the same field-model family.
2. A shape's editable state \(\mathbf{c}\) is separate from the field's learned weights \(\theta\).
3. Changing a named parameter changes \(\mathbf{c}\), not millions of weights and not a hand-authored alternate geometry backend.
4. Training or fitting maps the shape and its controls into this common representation.
5. Explicit equations, where available, generate trusted reference samples and test values. They are **oracles for supervision and validation**, not the production geometry representation.

For example, the sphere's radius is an input control named r. The field must respond to r through the same model used for the slab, cylinder, ellipsoid, SUBOFF, and later irregular shapes. We must not quietly switch the sphere back to an analytic signed-distance function at query time.

## 3. Field target and semantics

For the first prototype, use a truncated signed-distance target (TSDF), with a documented sign convention (negative inside, positive outside) and truncation distance \(\tau\):

\[
s(\mathbf{x}) = \operatorname{clip}(d(\mathbf{x}),-\tau,\tau)
\]

The model predicts this scalar field. The zero level set approximates the boundary. Do not assume the predicted value is an exact Euclidean distance everywhere; only claim an SDF if gradient-norm and distance-error tests support that claim.

The model must expose a stable query API:

- evaluate(points, state) -> field_values
- gradient(points, state) -> spatial_gradients
- set_state(name, value) or an immutable equivalent returning a new state
- extract_surface(state, resolution, bounds) -> mesh
- validate(state, reference) -> metrics

The API is common; the implementation must not route by shape type to different field evaluators.

## 4. Editability is a first-class requirement

A continuous field alone is insufficient if every design change requires retraining from scratch. The prototype therefore separates:

- **Shared model parameters** \(\theta\): learned field representation;
- **Explicit state** \(\mathbf{c}\): meaningful, named engineering controls such as radius, length, beam, draft, and local deformation handles;
- **Control semantics**: tests that specify what changing each control is expected to do.

A parameter edit should normally evaluate a new state through the same model without updating \(\theta\). For unfamiliar shapes, the research task is to learn or infer compact controls whose effects are predictable and local. If the prototype can only reproduce a shape but cannot edit it reliably, it has failed the central requirement.

Do not promise that every arbitrary topology change can be expressed by a small control vector. Record where the control model breaks down, and measure whether new shapes require fitting, a new state code, or a change to the shared model.

## 5. Staged experiment

### Stage A — controlled primitives

Use slab/plate, cylinder, ellipsoid, and sphere. Fit one shared conditional field model to generated samples across ranges of their state parameters. All queries go through the learned field, including shapes for which exact equations are known.

Tests:
- interpolation to parameter values not seen during training;
- direct radius/length/axis changes without weight updates;
- field and surface error against exact reference values;
- gradient direction and boundary-normal error;
- continuity and mesh convergence under resolution refinement;
- deterministic repeatability and inference cost.

A particularly important test: train over a range of sphere radii, then change r to an unseen value and query the same field model without retraining.

### Stage B — SUBOFF

Generate reference samples from the existing SUBOFF profile/geometry implementation, but do not use that implementation as the runtime evaluator. Fit the common field representation and expose meaningful controls (for example, length and radial scale, if the source geometry supports those semantics).

Compare field predictions against reference geometry and test whether parameter edits preserve the intended hull changes. The analytic/profile model remains a reference oracle only.

### Stage C — geometry without a convenient explicit formula

Use a fixed shared field architecture and the same query/extraction API. Fit from available mesh, point-cloud, or simulation-derived samples. Investigate compact controls and local deformation handles as **state inputs to the same field**, not as alternate geometry representations.

Evaluate whether edits can be made without retraining the shared weights, and whether locality, smoothness, and engineering constraints remain controllable. If they cannot, report the limitation rather than adding a second backend and calling it unified.

## 6. Required baselines and measurements

The baselines are for scientific comparison only; they are not runtime dispatch choices.

Report at minimum:

- boundary error (mean, 95th percentile, maximum);
- inside/outside classification accuracy;
- gradient/normal angular error;
- continuity and spurious-surface/self-intersection observations;
- edit error after changing each named control;
- locality: unintended displacement away from the edited region;
- whether a change required model retraining;
- parameter count, model size, inference throughput, training time, and surface extraction time;
- solver-facing mesh quality and repeatability.

Compare against known exact values for primitives and the available SUBOFF reference. These comparisons establish whether the unified representation is good enough; they do not justify retaining formula-based production paths.

## 7. First implementation boundary

Keep the existing lightweight TensorCAD core dependency-free. Put the research implementation behind an optional PyTorch extra and isolate it in a clearly named experimental module. Add a reproducible script that:

1. generates a fixed, seeded training/validation split over parameter ranges;
2. trains one conditional field model;
3. evaluates unseen states;
4. exports numerical metrics and extracted meshes;
5. runs editability tests with model weights frozen.

Do not claim success from attractive renderings alone. A successful first milestone must demonstrate measured accuracy **and** an unseen parameter edit without retraining.

## 8. Decision rule

Proceed with this unified field approach only if one model family can represent all Stage A shapes and SUBOFF with acceptable error, and meaningful controls work without weight updates for the tested edits. If it fails, diagnose the failure and revise the common representation or control parameterization. Do not respond to failure by silently dispatching to analytic formulas for some shapes and another representation for others.

The core research question is:

> Can a single continuous field, conditioned on a compact editable state, become the canonical geometry representation from regular primitives to complex engineering shapes, while explicit geometry serves as the test oracle?
