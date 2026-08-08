# Claim bundle — the Biological Claim Layer source

This directory is the **claim-layer source artifact** the `chemotactic-recruitment`
investigation realizes. It is a backend-neutral, typed bundle of graphs describing
*what is asserted* and *how it is (one way) made executable* — independent of the
CPM engine. It is the vertical-slice example for the Biological Claim Layer program
(see the ecosystem memory / RFC line): a future claim compiler would **lower this
bundle into** the three member studies rather than have them hand-authored.

| File | Layer | Answers |
|---|---|---|
| `claim.yaml` | Biological Claim Graph | What is asserted? (claim, context, intervention, observable, alternatives, applicability, falsifiers) |
| `mechanism.yaml` | Conceptual Mechanism Template | What can happen? (roles + capabilities + contribution + controls — no backend) |
| `realization.yaml` | Realization Graph | How is it implemented in viva-cpm? (bindings, fidelity, inherited defaults, exclusions, **semantic gap**, equivalence level) |

## Claim → studies (what the compiler would generate)

```
claim:chemotactic-recruitment
  ├─ baseline realization      → study:recruitment-baseline      (composite recruitment_baseline)
  ├─ + block-response          → study:recruitment-inhibited     (composite recruitment_inhibited)
  └─ + remove-cue (control)    → study:recruitment-adversarial   (composite recruitment_adversarial)
```

## The honest bit

`realization.yaml` records the **semantic gap** explicitly: the CPM chemotaxis
`lambda` is a *phenomenological* response strength (a bias on the copy-attempt
energy), not a receptor-level model. The realization is labelled
`phenomenological-equivalent` — it reproduces the recruitment phenotype and its
joint cue+response dependence, but is not mechanistically equivalent to a receptor
model. Alternative realizations (PhysiCell, receptor-kinetics ODE) that would
satisfy the *same* claim are listed as future backends.
