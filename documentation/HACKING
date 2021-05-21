## Hacking

We recommend using Guix to hack on this pipeline.  To enter a
reproducible environment with a known-good version of Guix use this
(slow) command:

```sh
USE_GUIX_INFERIOR=t guix environment --pure -m manifest.scm --preserve=GUIX_LOCPATH
```

To use your current Guix channels instead of the fixed set of
channels, just omit the `USE_GUIX_INFERIOR` shell variable:

```sh
guix environment --pure -m manifest.scm --preserve=GUIX_LOCPATH
```

To fetch code that is common to all PiGx pipelines run this:

```sh
git submodule update --init
```

Inside the environment you can then perform the usual build steps:

```sh
./bootstrap.sh # to generate the "configure" script
./configure
make
make check
```

To run the whole pipeline on the included test data sets execute this:

```sh
rm -r tests/output
make integration
```

You can run the pipeline from the source directory without installing
it by setting the `PIGX_UNINSTALLED` shell variable to any value.  This command runs the pipeline on the included test data:

```sh
PIGX_UNINSTALLED=t ./pigx-sars-cov2-ww -s tests/settings.yaml tests/sample_sheet.csv
```

Please, also see https://github.com/BIMSBbioinfo/pigx_sarscov2_ww/blob/main/documentation/user_installation_doc.md for setting up the pipeline and the needed databases.


## GitHub how-tos

in general:
- make small, atomic commits
- merge often to main branch
- make use of `git rebase` (in different situations, see below)


### if everything is set up
```
git add <file>
git commit -m 'message'
git push origin <branch name>
# consider making pull request
```

### contribute to someone else's repo assuming you have permission to push
```
git clone <link>
git checkout master
git checkout -b <new branch name e.g. fix/issue71>
# --> make some changes
git add <file>
git commit -m 'message'
git push origin <branch name>
# --> make pull request
```

### collaboration & rebase
in general:
- the one who contributes sth. should deal with merge conflicts before making a pull request (instead of the person merging the feature)
- resolving conflicts during a merge --> not easy, because one may not be the best person to decide on a resolution
- resolving conflicts during a rebase beforehand is easier 
- keep the commit history clean --> more useful

basics:
go to feature branch; fetch remote branches and then rebase local commits on top of the remote `main` branch
check if there are untracked changes, if so run `git stash` first
run `git stash pop` when the rebase is complete to get back the untracked changes
```
git fetch origin
git rebase origin/main
```
examplatory situation -- collaborating on one branch:
- Person A and B are working on the same branch X. Person A pushes her recent commits to X. Now, B has to git fetch origin + git rebase origin/X 
- if there are any, resolve conflicts (see https://codeinthehole.com/guides/resolving-conflicts-during-a-git-rebase/)

further reading on git workflows:
- a long version; note, that working with fork is not necessary when having push permission: https://medium.com/singlestone/a-git-workflow-using-rebase-1b1210de83e5
- a shorter version: https://gist.github.com/alper/846c187892fe2f74b78a


### modify commit history (example: squash)
```
git rebase -i HEAD~4 # HEAD~4 is an example and shows last 4 commits; -i for interactive
```
- git opens editor
- e.g. put 'squash' instead of 'pick' in front of the last commit, rearrange the order s.t. the last commit gets combined with the commit you want it to
- git again opens editor
- modify commit message
- if commit was already pushed before, do a force push to 'override' it
```
git push origin <branch name> --force
```


### take & apply commit from different branch
If one need just one commit from a different branch to continue developing, cherry-picking it might be useful
```
git cherry-pick <commit hash>
```
