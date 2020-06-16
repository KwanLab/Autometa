Contributing Guidelines
=======================

First off, thanks for taking the time to contribute! :clap::+1::tada:

#### Table of Contents

[Contributing Code](#contributing-code)

[How Can I Contribute?](#how-can-i-contribute)

* [Reporting Bugs](#reporting-bugs)
* [Suggesting Enhancements](#suggesting-enhancements)
* [Submitting pull requests](#pull-requests)

[Style guides](#style-guides)

* [git commit messages](#git-commit-messages)
* [documentation style guide](#documentation-style-guide)
* [python style guide](#python-style-guide)

### Contributing Code

1. Fork the KwanLab Autometa repository. You can do this [here](https://github.com/KwanLab/Autometa).
2. <details><summary>Add the KwanLab as the upstream remote.</summary><code>
    git remote add upstream https://github.com/KwanLab/Autometa.git</code></details>
3. Update your local repository with the most recent updates from the Kwan Lab.
  - <details><summary>Update your local repository.</summary><code>git pull upstream master</code></details>
  - <details><summary>Push changes to update your remote repository.</summary><code>git push origin master</code></details>
  - <details><summary>Ensure pre-commit hooks are installed.</summary><code>conda install pre-commit -y; pre-commit install</code></details>
4. <details><summary>Check out a branch corresponding to the feature you wish to add.</summary><code>git checkout -b your-new-feature master</code></details>
5. File an issue with the feature you plan on adding. This will open a channel of
communication with the core developers. See [suggesting enhancements](#suggesting-enhancements) for details.
6. Submit a pull request! To submit a pull request - see [instructions below](#pull-requests).


## How can I contribute?

### Reporting Bugs

> **Note:** If you find a **Closed** issue that seems like it is the same thing
that you're experiencing, open a new issue and include a link to the original
issue in the body of your new one.

When submitting a bug report, please follow the bug report [template](https://github.com/KwanLab/Autometa/blob/master/.github/ISSUE_TEMPLATE/bug_report.md).


### Suggesting Enhancements

Feature requests may be submitted by creating a new issue. When creating this issue, please follow [this template](https://github.com/KwanLab/Autometa/blob/master/.github/ISSUE_TEMPLATE/feature_request.md). After your feature request has been submitted, a maintainer will respond with a determination of whether this feature is appropriate for Autometa. If a maintainer has not responded within a reasonable time period, you can notify the Autometa team. See [notifying the team](#notifying-the-team) for details.

### Pull Requests

Pull requests have several goals:

- Maintain Autometa's quality
- Add features that are important to users
- Enable a sustainable system for Autometa's maintainers to review contributions

Please follow these steps to have your contribution considered by the maintainers:

1. Follow all instructions in respective [bug_fix](https://github.com/KwanLab/Autometa/blob/master/.github/PULL_REQUEST_TEMPLATE/bug_fix.md) or [feature_change](https://github.com/KwanLab/Autometa/blob/master/.github/PULL_REQUEST_TEMPLATE/feature_change.md) templates.
2. Follow the [style guides](#style-guides).
>Note: If writing a new file, a `template.py` file is provided within the autometa code base to help follow the structure requested by the Autometa team.
Please copy and rename this file before you start writing your feature.
The template file may be found [here](https://github.com/KwanLab/Autometa/blob/dev/docs/template.py), or you may find it within your cloned repository under 'Autometa/docs/template.py'

3. After you submit your pull request, verify that all [status checks](https://help.github.com/articles/about-status-checks/) are passing. <details><summary>What if the status checks are failing?</summary>If a status check is failing, and you believe that the failure is unrelated to your change, please leave a comment on the pull request explaining why you believe the failure is unrelated. A maintainer will re-run the status check for you. If we conclude that the failure was a false positive, then we will open an issue to track that problem with our status check suite.</details>

While the above must be satisfied, additional prerequisites may be present, depending on the type of pull request being issued (feature change or bug fix). The reviewer(s) may also ask you to complete additional changes before your pull request can finally be accepted.

## Style Guides

### Git Commit Messages

* Use the present tense ("Add feature" not "Added feature")
* Limit the first line to 72 characters or less
* Reference issues and pull requests liberally after the first line
* Consider starting the commit message with an applicable emoji:
  * :art: `:art:` when improving the format/structure of the code
  * :racehorse: `:racehorse:` when improving performance
  * :memo: `:memo:` when writing docs
  * :penguin: `:penguin:` when fixing something on Linux
  * :apple: `:apple:` when fixing something on macOS
  * :bug: `:bug:` when fixing a bug
  * :fire: `:fire:` when removing code or files
  * :green_heart: `:green_heart:` when fixing the CI build
  * :white_check_mark: `:white_check_mark:` when adding tests
  * :arrow_up: `:arrow_up:` when upgrading dependencies
  * :arrow_down: `:arrow_down:` when downgrading dependencies

### Documentation Style Guide

Documentation is hosted on reathedocs.org which uses sphinx, a python documentation generator. Therefore a common syntax for
documenting files within the source code is required. The Autometa team follows numpy syntax to automate as much of the documentation build as possible.

* Numpy documentation [style guide](https://numpydoc.readthedocs.io/en/latest/format.html)
* Sphinx restructured text [style guide](http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html#paragraphs)

### Python Style Guide

We have moved all development to using the [black](black.readthedocs.io) formatter and run
[pre-commit](https://pre-commit.com/) hooks with black such that all of the fuss with formatting according to code
specifications can be ignored! When first cloning the repository, you should first install
pre-commit and run `pre-commit install`. This will ensure any commits you make will be
formatted appropriately.

### Notifying the Team

You may notify the Autometa team corresponding to the respective update/bug by mentioning them in a comment within your issue or pull request.

>Note: Please only use <code>@KwanLab/autometa</code> if you receive no response from any of the other teams below.

| Team                                   | Area of Development                       |
| -------------------------------------- | ----------------------------------------- |
| <code>@KwanLab/autometa-core</code>    | Members developing the core functionality |
| <code>@KwanLab/autometa-website</code> | Members developing the website            |
| <code>@KwanLab/autometa</code>         | All Members of the Autometa team          |
