import re

import pytest
from playwright.sync_api import Locator, Page, expect


def authentication_scenarios(*args):
    if len(args) == 2:
        # Run only the authenticated scenario if "authenticated" mark is used, otherwise run both scenarios
        return [
            pytest.param(True, marks=[pytest.mark.authenticated]),
            pytest.param(False),
        ]
    else:
        # Run the test regardless of the "authenticated" mark
        return [pytest.param(args[0], marks=[pytest.mark.authenticated])]


def ensure_accepted_cookies(page: Page):
    if page.get_by_text(
        "This website uses cookies to ensure you get the best experience on our website. "
    ).is_visible():
        with page.expect_popup() as privacy_popup_info:
            page.get_by_text("Learn more").click()
        privacy_policy_page = privacy_popup_info.value

        expect(privacy_policy_page).to_have_url(
            "https://www.genialis.com/privacy-policy/"
        )
        privacy_policy_page.close()

        page.get_by_label("Accept cookies").click()

    expect(
        page.get_by_text(
            "This website uses cookies to ensure you get the best experience on our website. "
        )
    ).not_to_be_visible()


def log_in(page: Page, user: dict):
    page.wait_for_load_state("domcontentloaded")
    expect(page.get_by_label("Fetching in progress")).not_to_be_visible()
    page.wait_for_load_state("domcontentloaded")
    if page.get_by_role("button", name="Login").is_visible():
        page.get_by_role("button", name="Login").click()

    if page.get_by_text(
        "You're being redirected to a SSO login page. Please click the button below"
    ).is_visible():
        page.get_by_role("button", name="Log in").click()

    page.get_by_label("Email address").fill(user["username"])
    page.get_by_label("Password").fill(user["password"])
    page.get_by_role("button", name="Continue", exact=True).click()
    expect(page.get_by_label("Fetching in progress")).not_to_be_visible()


def add_genes_to_selection(page: Page, genes):
    for gene in genes:
        page.get_by_placeholder("Search for a gene").fill(gene)
        page.get_by_role("option", name=f"{gene}").click()

    for gene in genes:
        expect(page.get_by_role("button", name=gene)).to_be_visible()

    page.get_by_text("GenesHistory").click()


def wait_for_widgets_loaded(page: Page, widgets):
    if "Gene Ontology Enrichment" in widgets:
        expect(page.get_by_text("Enriched terms not found.")).not_to_be_visible(
            timeout=60000
        )

    for widget in widgets:
        expect(
            page.locator("div")
            .filter(has_text=re.compile(rf"^{widget}$"))
            .get_by_role("progressbar")
        ).not_to_be_visible(timeout=60000)
        expect(
            page.locator("div")
            .filter(has_text=re.compile(rf"^{widget}$"))
            .get_by_test_id("ScheduleIcon")
        ).not_to_be_visible(timeout=60000)


def full_layout_locator(page: Page) -> Locator:
    return (
        page.locator("div")
        .filter(
            has_text="Time series and Gene SelectionFilter time seriesFilter time series Project 1 Det"
        )
        .nth(1)
    )


def widget_locator(page: Page, widget_name: str) -> Locator:
    return page.locator("div.react-grid-item").filter(
        has=page.locator("div").filter(has_text=re.compile(rf"^{widget_name}$"))
    )
