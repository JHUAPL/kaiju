Contributing Guide
==================

Introduction
------------

This guide assumes that you have already checked out a local copy of our git
repository using the other pages in our wiki. If not, please check other pages
first.

This guide will talk about concepts such as pull requests and branching.
Details about them are beyond the scope of this guide, so for additional
information about them please refer to resources such as these:

* `Pull Request Tutorial <https://support.atlassian.com/bitbucket-cloud/docs/tutorial-learn-about-bitbucket-pull-requests/>`_
* `Making A Pull Request <https://www.atlassian.com/git/tutorials/making-a-pull-request>`_
* `Git Branching Overview <https://www.atlassian.com/git/tutorials/learn-branching-with-bitbucket-cloud>`_

Branch Setup
------------

In our repository we have designated the branch named "master" as both the
"Main" branch, and as the "Production" branch.

Because it is the "Main" branch, whenever anyone downloads our repository,
this is the default branch that they will see and use.

Because it is the "Production" branch, we only merge new changes and features
into it periodically once we are confident that these new features are stable
and correct.

We have also designated the branch named "development" as the "Development"
branch. This branch is where new and experimental features should be added so
that they can be tested thoroughly.

Periodically, the "development" branch gets pulled into the "master" branch
and marks a new release of the code. This way anyone using the code does not
accidentally check out an unstable or untested feature (but if they want to,
they can use the "development" branch and all of its latest features).

Branch Access
-------------

The "development" and "master" branches are locked from direct access, so
noone is allowed to push any changes directly to either branch. All changes
are required to be done by pull requests (explained below). This ensures that
all code added to the "development" branch has some level of review, and
prevents anyone from accidentally pushing changes directly to the production
"master" branch.

Contributing Code
-----------------

Whenever someone has a new feature, bugfix, or anything that they want to add
to the repository, they begin by making a new branch off of the "development"
branch. We recommend naming the new branch something descriptive about the
work being added, but it can be named whatever you like. This new branch
**must** come off of the "development" branch, because that is where new work
gets merged into the code base.

Once the new feature has been completed (and **tested**\ ), you can submit a
pull request through the bitbucket website so that this feature branch can be
merged into the "development" branch. One or more developers or other
appropriate experts will review the changes in the branch, and then either
approve the pull request (at which point the branch is merged into
"development" and the code is now part of the repository), or reject the pull
request with an explanation as to what additional changes may be needed to the
code in the feature branch.

A list of recommended reviewers is available
:doc:`here <recommendedReviewers>`.
