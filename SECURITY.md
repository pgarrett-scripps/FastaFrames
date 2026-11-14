# Security policy

## Supported versions

Security fixes go into the latest release of fastaframes on PyPI.

## Reporting a vulnerability

Please report vulnerabilities privately through GitHub security advisories:
open the repository's **Security** tab and choose **Report a vulnerability**
(<https://github.com/pgarrett-scripps/FastaFrames/security/advisories/new>).
Do not open a public issue.

Include the affected version, steps to reproduce, and the impact you expect.
You should get a reply within seven days. Once a fix is released, the advisory
is published and reporters are credited if they wish.

FastaFrames parses files that come from outside your control. A FASTA file that
makes it raise an unexpected error, consume unbounded memory, or write data it
was not given is in scope. Note that nearly all parsing happens in
[fastatacular](https://github.com/tacular-omics/fastatacular); if the problem is
in the parser rather than in the DataFrame mapping, we will forward the report
there and say so in the advisory.
