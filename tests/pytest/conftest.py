"""
Custom pytest configuration for DARTassembler tests.

Adds two custom categories:
- SMALLCHANGE: numerically different but chemically identical (counts as success).
- BIGCHANGE: large or relevant changes (shown as failures in local tools like PyCharm,
  but the GitHub CI will still return exit code 0 if BIGCHANGE is the only failure type).
"""
import pytest

def pytest_addoption(parser):
    grp = parser.getgroup("change categories")
    grp.addoption(
        "--bigchange-policy",
        action="store",
        default="fail",
        choices=("allow", "fail", "warn"),
        help=(
            "How to treat tests labelled BIGCHANGE via request.node.bigchange:\n"
            "  allow -> force pass\n"
            "  fail  -> force fail (even if assertions passed)\n"
            "  warn  -> keep real outcome, only relabel"
        ),
    )

def pytest_configure(config):
    # Nothing to register as a marker anymore, but we expose a helper fixture.
    pass

# Optional ergonomic helper so you don't have to touch request.node directly
@pytest.fixture
def change(request):
    class _ChangeSetter:
        def __init__(self, node):
            self._node = node
        def big(self, reason: str = ""):
            self._node.bigchange = True
            self._node.bigchange_reason = reason
            # Clear the opposite if user mixed calls
            if hasattr(self._node, "smallchange"):
                delattr(self._node, "smallchange")
        def small(self, reason: str = ""):
            self._node.smallchange = True
            self._node.smallchange_reason = reason
            if hasattr(self._node, "bigchange"):
                delattr(self._node, "bigchange")
        def set(self, kind: str):
            kind = kind.lower()
            if kind.startswith("big"):
                self.big()
            elif kind.startswith("small"):
                self.small()
            else:
                raise ValueError("kind must be 'big' or 'small'")
    return _ChangeSetter(request.node)

@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    # run test and get the original report
    outcome = yield
    rep = outcome.get_result()

    # only touch the test body (not setup/teardown)
    if rep.when != "call":
        return

    is_big = bool(getattr(item, "bigchange", False))
    is_small = bool(getattr(item, "smallchange", False))

    # If both were set by mistake, treat as BIGCHANGE (more conservative)
    if is_big and is_small:
        is_small = False
        rep.longrepr = rep.longrepr or "Both bigchange and smallchange set; treating as BIGCHANGE."

    policy = item.config.getoption("--bigchange-policy")

    # SMALLCHANGE: always pass, custom label
    if is_small:
        rep.smallchange = True
        rep.outcome = "passed"
        rep.longrepr = None  # suppress any prior failure output
        rep.change_reason = getattr(item, "smallchange_reason", "Minor change")
        return

    # BIGCHANGE: policy-driven behavior
    if is_big:
        rep.bigchange = True
        rep.change_reason = getattr(item, "bigchange_reason", None)
        if policy == "allow":
            rep.outcome = "passed"
            rep.longrepr = None
        elif policy == "fail":
            rep.outcome = "failed"
            if not rep.longrepr:
                rep.longrepr = "BIGCHANGE not allowed (use --bigchange-policy=allow to accept)"
        elif policy == "warn":
            # keep real outcome, only relabel in pytest_report_teststatus
            pass

def pytest_report_teststatus(report, config):
    # Give custom progress letter + verbose label for our categories (call phase only)
    if report.when != "call":
        return

    if getattr(report, "smallchange", False):
        # Count as passed; show S / SMALLCHANGE
        reason = getattr(report, "change_reason", "")
        verbose = f"SMALLCHANGE – {reason}" if reason else "SMALLCHANGE"
        return ("passed", "S", verbose)

    if getattr(report, "bigchange", False):
        policy = config.getoption("--bigchange-policy")
        reason = getattr(report, "change_reason", "")
        verbose = f"BIGCHANGE[{policy}] – {reason}" if reason else f"BIGCHANGE[{policy}]"
        # Keep summary category aligned with true outcome for exit codes
        return (report.outcome, "B", verbose)