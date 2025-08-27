import os

import pytest
from _pytest.fixtures import SubRequest
from playwright.sync_api import Browser, BrowserContext, Page, Playwright, expect
from pytest_metadata.plugin import metadata_key
from tests.utils import add_genes_to_selection, ensure_accepted_cookies, log_in
from urls import get_urls

EXPRESSIONS_E2E_USERNAME = "EXPRESSIONS_E2E_USERNAME"
EXPRESSIONS_E2E_PASSWORD = "EXPRESSIONS_E2E_PASSWORD"


# If variable with this name is set to True traces are stored.
CI = "CI"


# Set test variables

page_with_context: Page = None


def _get_env_var(varname: str) -> str:
    value = os.getenv(varname)
    assert value, f"{varname} is not set"
    return value


# set up environment


def pytest_addoption(parser: pytest.Parser):
    parser.addoption("--env")


@pytest.fixture(scope="session", autouse=True)
def env(request: pytest.FixtureRequest):
    return request.config.getoption("--env")


def pytest_configure(config: pytest.Config):
    config.addinivalue_line(
        "markers", "authenticated: Runs tests only for authenticated users."
    )


@pytest.fixture(scope="session", autouse=True)
def urls(env):
    return get_urls(env)


@pytest.fixture(scope="session", autouse=True)
def user():
    return {
        "username": _get_env_var(EXPRESSIONS_E2E_USERNAME),
        "password": _get_env_var(EXPRESSIONS_E2E_PASSWORD),
    }


# Set up context


@pytest.fixture(scope="session", autouse=True)
def set_up_session(
    playwright: Playwright,
    user: dict,
    urls,
    request: pytest.FixtureRequest,
):
    print("Setup started.")

    browser = playwright.chromium.launch(headless=True, slow_mo=500)
    context = browser.new_context(ignore_https_errors=True)
    context.tracing.start(screenshots=True, snapshots=True, sources=True)
    page = context.new_page()

    try:
        page.goto(urls["root"])
        request.config.stash[metadata_key]["Base URL"] = urls["base"]

        ensure_user_logged_in(page, user)
        ensure_accepted_cookies(page)

        page.context.storage_state(path="state.json")
    except Exception as e:
        context.tracing.stop(path="artifacts/traces/setup_session.zip")
        raise e

    yield

    context.close()
    browser.close()


@pytest.fixture()
def new_context(authenticated, name_of_test, browser: Browser):
    try:
        tracing = _get_env_var(CI) == "true"
    except AssertionError:
        tracing = False
    context = browser.new_context(
        storage_state=("state.json" if authenticated else None)
    )
    context.grant_permissions(["clipboard-read", "clipboard-write", "notifications"])

    if tracing:
        context.tracing.start(screenshots=True, snapshots=True, sources=True)

    global page_with_context
    page_with_context = context.new_page()

    yield context

    if tracing:
        context.tracing.stop(path=f"artifacts/traces/{name_of_test}.zip")

    context.close()


@pytest.fixture()
def authenticated(request: SubRequest):
    return request.param


@pytest.fixture()
def name_of_test(request: SubRequest):
    """Get the name of the test."""
    return request.node.name


@pytest.fixture()
def root_page(new_context: BrowserContext, urls):
    page: Page = new_context.new_page()
    page.goto(urls["root"])
    ensure_accepted_cookies(page)

    yield page


@pytest.fixture()
def root_page_with_gene_selection(new_context: BrowserContext, urls, genes):
    page: Page = new_context.new_page()
    page.goto(urls["root"])
    ensure_accepted_cookies(page)
    page.wait_for_timeout(500)
    add_genes_to_selection(page, genes)

    yield page


@pytest.fixture()
def widgets():
    return [
        "Time series and Gene Selection",
        "Expression Time Courses",
        "Differential expressions",
        "Hierarchical Clustering",
        "Gene Ontology Enrichment",
    ]


@pytest.fixture()
def genes():
    return ["aplA", "gpaA", "tifA", "enoB"]


def ensure_user_logged_in(page: Page, user: dict):
    log_in(page, user)
    expect(page.get_by_role("button", name="Login")).not_to_be_visible()
