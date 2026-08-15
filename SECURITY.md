# Security Policy

## Reporting a Vulnerability

Please report security vulnerabilities by opening a new [GitHub security
advisory](https://github.com/MDAnalysis/membrane-curvature/security/advisories/new).

> [!WARNING]
> Please do **not** file a public GitHub issue for security reports.

When reporting, if possible, include:
- A clear description of the issue and potential impact.
- Steps to reproduce or a proof of concept.
- Affected MembraneCurvature version and environment details.
- Any suggested mitigations.

You can also send an email to `mdanalysis@numfocus.org`, which is an alias to
a subset of the MDAnalysis maintainers' team.

If the security vulnerability is accepted, a patch will be crafted privately
in order to prepare a dedicated bugfix release as timely as possible (depending
on the complexity of the fix).

## Supported Versions

Only the main branch (the 2.x release line) receives updates, including
security fixes. Older release lines (e.g., 1.1.x) are not maintained. If you are
using an older version, please upgrade to the latest 2.x release to receive
fixes.

When reporting a potential vulnerability, confirm that it reproduces against
the latest 2.x version.

| Version       | Supported          |
| ------------- | ------------------ |
| 2.0.0         | ✅ |
| < 1.1.2       | ❌ |
