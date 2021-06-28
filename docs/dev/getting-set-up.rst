Getting set up
==============

Working with `git`
------------------

This section of the docs covers our practices for working with `git` on our codebase. For more in-depth guides, we can recommend a few sources:

For a more complete git tutorials we recommend checking out:

* `Atlassian's git tutorial <https://www.atlassian.com/git/tutorials>`__ -- Beginner friendly introductions to the git command line interface
* `Setting up git for GitHub <https://docs.github.com/en/free-pro-team@latest/github/getting-started-with-github/set-up-git>`__ -- Configuring git to work with your GitHub user account

.. _forking-and-cloning:

Forking and cloning
~~~~~~~~~~~~~~~~~~~

To get the code, and be able to push changes back to the main project, you'll need to (1) fork the repository on github and (2) clone the repository to your local machine.

This is very straight forward if you're using `GitHub's CLI <https://cli.github.com>`__:

.. code:: shell

    $ gh repo fork theislab/scanpy --clone --remote

This will fork the repo to your github account, create a clone of the repo on your current machine, add our repository as a remote, and set the `master` development branch to track our repository.

To do this manually, first make a fork of the repository by clicking the "fork" button on our main github package. Then, on your machine, run:

.. code:: shell

    # Clone your fork of the repository (substitute in your username)
    git clone https://github.com/{your-username}/scanpy.git
    # Enter the cloned repository
    cd scanpy
    # Add our repository as a remote
    git remote add upstream https://github.com/theislab/scanpy.git
    # git branch --set-upstream-to "upstream/master"

`pre-commit`
~~~~~~~~~~~~

We use `precommit <https://pre-commit.com>`__ to run some styling checks in an automated way.
We also test against these checks, so make sure you follow them!

You can install pre-commit with:

.. code:: shell

    pip install pre-commit

You can then install it to run while developing here with:

.. code:: shell

    pre-commit install

From the root of the repo.

If you choose not to run the hooks on each commit, you can run them manually with `pre-commit run --files={your files}`.

.. _creating-a-branch:

Creating a branch for your feature
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All development should occur in branches dedicated to the particular work being done.
Additionally, unless you are a maintainer, all changes should be directed at the `master` branch.
You can create a branch with:

.. code:: shell

    git checkout master                 # Starting from the master branch
    git pull                            # Syncing with the repo
    git checkout -b {your-branch-name}  # Making and changing to the new branch

.. _open-a-pr:

Open a pull request
~~~~~~~~~~~~~~~~~~~

When you're ready to have your code reviewed, push your changes up to your fork:

.. code:: shell

    # The first time you push the branch, you'll need to tell git where
    git push --set-upstream origin {your-branch-name}
    # After that, just use
    git push

And open a pull request by going to the main repo and clicking *New pull request*.
GitHub is also pretty good about prompting you to open PRs for recently pushed branches.

We'll try and get back to you soon!

.. _dev-environments:

Development environments
------------------------

It's recommended to do development work in an isolated environment.
There are number of ways to do this, including conda environments, virtual environments, and virtual machines.

We think the easiest is probably conda environments. Simply create a new environment with a supported version of python and make a :ref:`development install <dev-install-instructions>` of `scanpy`.
