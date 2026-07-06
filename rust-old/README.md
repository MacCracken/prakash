# rust-old — archived Rust source (pre-2.0.0)

This directory holds the original **Rust** implementation of prakash (crate
version 1.2.0), preserved for reference while the library is ported to the
[Cyrius](https://github.com/MacCracken/cyrius) language.

It is a self-contained Cargo project: `Cargo.toml` + `src/` + `tests/` +
`benches/` + `examples/` + `scripts/` moved here verbatim. From this directory,
`cargo build` / `cargo test` still work with the pinned toolchain
(`rust-toolchain.toml`).

**It is the source of truth for the port** — each Cyrius module in the repo
root's `src/*.cyr` is translated from the corresponding file here, with the
ported Rust tests becoming `tests/*.tcyr` assertions. See
[`docs/development/cyrius-port-plan.md`](../docs/development/cyrius-port-plan.md).

Once the port reaches feature parity and ships as prakash **2.0.0**, this
directory is removed (the Rust history remains available in the pre-2.0 git
tags), mirroring how `hisab` retired its Rust source at its own 2.0.0 port.
