# How to contribute

We welcome contributions from external contributors, and this document
describes how to merge code changes into this membrane_curvature.

## Getting Started

* Make sure you have a [GitHub account](https://github.com/signup/free).
* [Fork](https://help.github.com/articles/fork-a-repo/) this repository on GitHub.
* On your local machine,
  [clone](https://help.github.com/articles/cloning-a-repository/) your fork of
  the repository.

## Development Environment

We strongly recommend using [uv](https://docs.astral.sh/uv/) as a development environment for MembraneCurvature. To set up a `uv` dev environment, run the following commands from the root of the repository:

```bash
uv sync --group dev
uv run pre-commit install
```

### Running tests

Tests are run with [pytest](https://docs.pytest.org/en/stable/) and are a fundamental part of the development process. All tests must pass before a PR can be merged.

To run the tests:

```bash
uv run pytest
```

To get the coverage report, run the tests with the `--cov` flag:

```bash
uv run pytest membrane_curvature/tests  --cov=membrane_curvature --cov-report=html --cov-report=term
```

### Building the documentation

To build the documentation, run the following command from the root of the repository:

```bash
uv run sphinx-build -b html docs/ build/html
```

For more information on building the documentation, see the [docs/README.md](../docs/README.md) file.

## Making Changes

* Add some really awesome code to your local fork.  It's usually a [good
  idea](http://blog.jasonmeridth.com/posts/do-not-issue-pull-requests-from-your-master-branch/)
  to make changes on a
  [branch](https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/)
  with the branch name relating to the feature you are going to add.
* Add the relevant changes to the [CHANGELOG.rst](../CHANGELOG.rst) file and add your name to the [AUTHORS](../AUTHORS) file if it's your first contribution.
* When you are ready for others to review and comment on your new feature,
  navigate to your fork of membrane_curvature on GitHub and open a [pull
  request](https://help.github.com/articles/using-pull-requests/) (PR). Note that
  after you launch a PR from one of your fork's branches, all
  subsequent commits to that branch will be added to the open pull request
  automatically.  Each commit added to the PR will be validated for
  mergability, compilation and test suite compliance; the results of these tests
  will be visible on the PR page.
* If you're providing a new feature, you must add test cases and documentation.
* When the code is ready to go, make sure you run the test suite using pytest.
* When you're ready to be considered for merging, check the "Ready to go"
  box on the PR page to let the membrane_curvature devs know that the changes are complete.
  The code will not be merged until this box is checked, the continuous
  integration returns checkmarks,
  and multiple core developers give "Approved" reviews.

> [!NOTE]
> All tests must pass before a PR can be merged. CI runs the test suite, linters, and type checks.
> New code warrants new tests, and all added code should be covered by the test suite.
> Coverage is reported to Codecov and it is reported in the PR. Maintainers may block merge if code coverage drops below 100%.

# Additional Resources

* [General GitHub documentation](https://help.github.com/)
* [PR best practices](http://codeinthehole.com/writing/pull-requests-and-other-good-practices-for-teams-using-github/)
* [A guide to contributing to software packages](http://www.contribution-guide.org)
* [Thinkful PR example](http://www.thinkful.com/learn/github-pull-request-tutorial/#Time-to-Submit-Your-First-PR)
