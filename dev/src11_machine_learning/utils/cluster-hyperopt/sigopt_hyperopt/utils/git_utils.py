import git


class GitException(Exception):
    """Exception thrown in the domain of git"""
    pass


def checkout_version(repo: git.repo.base.Repo, version: str, git_uri: str):
    try:
        repo.git.checkout(version)
    except git.exc.GitCommandError:
        raise GitException(f"Unable to checkout version {version} of repository {git_uri}. "
                           f"Please ensure that the version exists in the repo!")


def checkout_branch(repo: git.repo.base.Repo, origin: git.remote.RemoteReference, branch: str, git_uri: str):
    refs = repo.refs
    # refs are saved like origin/branch
    if branch not in [ref.name.split("/")[-1] for ref in refs]:
        raise GitException(f"The branch {branch} could not be found in repository "
                           f"{git_uri}. Please ensure that the branch exists in the repo!")
    repo.create_head(branch, origin.refs[branch])
    repo.heads[branch].checkout()
