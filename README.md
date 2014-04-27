E80-6DOF
========

6 degree of freedom model for rocket flight


New to Github? (original Author: Joshua Vasquez)
==============

Linnea Shin has compiled a [wonderful guide of tips for working with Github](https://gist.github.com/dunvi/5080965).  Not all of them
may be 100% correct, since the file is old and relates directly to Apple's XCode IDE, but it's
worth a read.


The following blurb below is a chunk of notes that I wrote while Sarah Lichtman and I were getting started
on a big Github project.
It covers a basic workflow with Git.  Git feels a little counterintuitive at first, so I did my best
to document normal usage.

## Basic Git Usage
1. Always pull from the remote repository first, to make sure you have the latest code.
1. For large changes, work on a distinct [branch](gitref.org/branching/).  When your done, and want to save changes, you can either 
[save that branch up on the remote repository](gitready.com/beginner/2009/02/02/push-and-delete-branches.html) or merge the working branch with the master and update the remote repository.

Basic workflow will look like:


Pull down the latest changes:
```
git pull
```

Make a new branch and switch to it with checkout:
```
git branch cool_new_features 
git checkout cool_new_features
```

Do some work on that branch...

When you're ready to commit, you have two options:

###Option 1: push up the new branch that you're working on:
```
git add .
git commit -m 'just started working on cool_new_features in a new branch'
git push origin cool_new_features
```
Git will save all of your branches in the remote repository if you add them.
Take care though: if a more recent version of cool_new_features is sitting in 
the remote repository, just like with the master branch, you'll have to git 
pull to merge the more recent additions into yours before you can push up 
your changes. Furthermore, since you're now dealing with multiple branches, 
you must specify where you're pulling from. Thus, a scenario like this would 
look like:

```
git pull origin cool_new_features
git add .
git commit -m 'just started working on cool_new_features in a new branch'
git push origin cool_new_features
```


###Option 2: merge with the master branch on the remote repository:
```
git checkout master
git pull
git merge cool_new_features 
git add .
git commit 'merged cool_new_features into master'
git push origin master
```
Now, you can go back to working on the same branch by checking it out, or you can delete it altogether. When  you're done working on a branch, you can delete it with:
```
git push origin --delete cool_new_features
```
Finally, to delete that branch locally on your machine, 
```
git branch -d cool_new_features
```
