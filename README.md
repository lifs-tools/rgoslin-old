# R implementation for parsing of lipid shorthand nomenclature names, version 2.0
[![Build Status](https://travis-ci.org/lifs-tools/rgoslin.svg?branch=master)](https://travis-ci.org/lifs-tools/rgoslin)[![codecov](https://codecov.io/gh/lifs-tools/rgoslin/branch/master/graph/badge.svg)](https://codecov.io/gh/lifs-tools/rgoslin)[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3757672.svg)](https://doi.org/10.5281/zenodo.3757672)

This project is a parser, validator and normalizer implementation for shorthand lipid nomenclatures, base on the Grammar of Succinct Lipid Nomenclatures project.

[https://github.com/lifs-tools/goslin](Goslin) defines multiple grammers compatible with ANTLRv4 for different sources of shorthand lipid nomenclature. This allows to generate parsers based on the defined grammars,
which provide immediate feedback whether a processed lipid shorthand notation string is compliant with a particular grammar, or not.

Here, rgoslin 2.0 uses the Goslin grammars and the cppgoslin parser to support the following general tasks:

1. Facilitate the parsing of shorthand lipid names dialects.
2. Provide a structural representation of the shorthand lipid after parsing.
3. Use the structural representation to generate normalized names.

## Related Projects

- [This project](https://github.com/lifs-tools/rgoslin)
- [Goslin grammars and reference test files](http://github.com/lifs-tools/goslin)
- [C++ implementation](https://github.com/lifs-tools/cppgoslin)
- [Java implementation](https://github.com/lifs-tools/jgoslin)
- [Python implementation](https://github.com/lifs-tools/pygoslin)
- [Webapplication and REST API](https://github.com/lifs-tools/goslin-webapp)

## Installation ##

### Prerequisites
Install the `devtools` package with the following command.
```R
if(!require(devtools)) { install.packages("devtools") }
```
  
Run

```R
  install_github("lifs-tools/rgoslin")
```
to install from the github repository.

This will install the latest, potentially unstable development version of the package with all required dependencies into your local R installation.

If you want to use a proper release version, referenced by a Git tag (here: `v2.0.1`) install the package as follows:

```R
  install_github("lifs-tools/rgoslin", ref="v2.0.1")
```

If you want to work off of a specific branch (here: `adding_masses`), install the package as follows:

```R
  install_github("lifs-tools/rgoslin", ref="adding_masses")
```

If you also want to build the help and vignette, add the following arguments:

```R
  install_github("lifs-tools/rgoslin", ref="adding_masses", build_manual = TRUE, build_vignettes = TRUE)
```

If you have cloned the code locally, use devtools as follows.
Make sure you set the working directory to where the API code is located.
Then execute

```R
library(devtools)
install(".")
```

To run the tests, execute
```R
library(devtools)
test()
```

### Usage

To load the package, start an R session and type

```R
  library(rgoslin)
```

Type the following to see the package vignette / tutorial:

```R
  vignette('introduction', package = 'rgoslin')
```

## Adding cppgoslin as a Git subtree

In the root of your git project, run the git subtree command, with `<PREFIX>` replaced by the subdirectory path where you want the subtree to live:

~~~~
git subtree add --prefix=<PREFIX> https://github.com/lifs-tools/cppgoslin.git master
~~~~

Instead of `master`, you can choose any other branch or tag to clone.
For more information on git subtree, see [Git Subtree](https://github.com/git/git/blob/master/contrib/subtree/git-subtree.txt) or [this article](https://blog.developer.atlassian.com/the-power-of-git-subtree/).


## Pulling and pushing of a Git subtree
For pulling and pushing, you have to change into the root directory of the host repository and execute the following commands:

### Pulling
~~~~
git subtree pull --prefix=<PREFIX> https://github.com/lifs-tools/cppgoslin.git master
~~~~

### Pushing
~~~~
git subtree push --prefix=<PREFIX> https://github.com/lifs-tools/cppgoslin.git master
~~~~

Alternatively, you can create shortcuts/aliases in your repository's `.git/config` file:

~~~~
[alias]
    # the acronym stands for "subtree pull"
    cppgoslin-pull = "!f() { git subtree pull --prefix <PREFIX> git@github.com:lifs-tools/cppgoslin.git master; }; f"
    # the acronym stands for "subtree push"
    cppgoslin-push = "!f() { git subtree push --prefix <PREFIX> git@github.com:lifs-tools/cppgoslin.git master; }; f"
~~~~

Make sure to replace `<PREFIX>` with the proper path from your repository root directory to the directory where you placed your subtree in!

This allows you to run `git cppgoslin-pull` to pull the latest master version, or `git cppgoslin-push` to push your latest local commits on the cppgoslin subtree to the upstream repository.

