from datetime import datetime
from pathlib import Path

import pytest
from playwright.sync_api import Page, expect
from tests.utils import (
    authentication_scenarios,
    wait_for_widgets_loaded,
    widget_locator,
)


@pytest.mark.parametrize(
    "authenticated", authentication_scenarios(True, False), indirect=True
)
def test_landing_page(root_page, urls, widgets):
    page: Page = root_page

    page.get_by_role("link", name="dictyExpress logo").click()
    expect(page).to_have_url(urls["base"])
    expect(
        page.get_by_text(
            "An interactive, exploratory data analytics web-app that provides access to gene "
        )
    ).to_be_visible()

    page.get_by_role("link", name="Run dictyExpress").click()

    expect(page).to_have_url(urls["root"])

    for widget in widgets:
        expect(page.get_by_text(widget)).to_be_visible()


@pytest.mark.parametrize(
    "authenticated", authentication_scenarios(False), indirect=True
)
def test_layout(root_page_with_gene_selection, widgets, assert_snapshot):
    page: Page = root_page_with_gene_selection

    wait_for_widgets_loaded(page, widgets)

    page.get_by_text("Time series and Gene Selection").click()

    original_widget_styles = [
        widget_locator(page, widget).get_attribute("style") for widget in widgets
    ]

    for widget in widgets:
        widget_resize_btn = (
            page.locator(".react-resizable")
            .filter(has_text=widget)
            .locator(".react-resizable-handle")
        )
        btn_box = widget_resize_btn.bounding_box()

        widget_resize_btn.hover()
        page.mouse.down()
        page.mouse.move(btn_box["x"] - 50, btn_box["y"] - 50)
        page.wait_for_timeout(500)
        page.mouse.up()

    changed_widget_styles = [
        widget_locator(page, widget).get_attribute("style") for widget in widgets
    ]

    assert original_widget_styles != changed_widget_styles

    page.get_by_role("button", name="Default layout").click()
    page.wait_for_timeout(1000)

    default_widget_styles = [
        widget_locator(page, widget).get_attribute("style") for widget in widgets
    ]

    assert original_widget_styles == default_widget_styles


@pytest.mark.parametrize(
    "authenticated", authentication_scenarios(True, False), indirect=True
)
def test_bookmark_link(root_page_with_gene_selection, genes, widgets, assert_snapshot):
    page: Page = root_page_with_gene_selection

    wait_for_widgets_loaded(page, widgets)

    page.get_by_role("button", name="Bookmark").click()
    expect(page.get_by_role("heading", name="Bookmark URL")).to_be_visible()
    url = page.get_by_role("link").get_attribute("href")

    page.get_by_role("button", name="Copy").click()
    clipboard_content = page.evaluate("navigator.clipboard.readText()")
    assert clipboard_content == url

    with page.expect_popup() as bookmark_info:
        page.get_by_role("link", name=url).click()
    bookmark_page = bookmark_info.value

    expect(bookmark_page).to_have_url(url)

    expect(
        bookmark_page.get_by_text("Computing Gene Ontology Enrichment for 4 genes.")
    ).to_be_visible(timeout=10000)

    bookmark_page.get_by_text("GenesHistory").click()
    wait_for_widgets_loaded(bookmark_page, widgets)

    for gene in genes:
        expect(bookmark_page.get_by_role("button", name=gene)).to_be_visible()

    expect(bookmark_page.get_by_label("Differential expression")).to_contain_text(
        "D. discoideum (prespore vs. prestalk)"
    )
    expect(bookmark_page.get_by_label("Clustering Linkage")).to_contain_text("Average")
    expect(bookmark_page.get_by_label("Distance Measure")).to_contain_text("Euclidean")
    expect(bookmark_page.get_by_label("Aspect")).to_contain_text("Biological process")
    expect(bookmark_page.get_by_label("p-value")).to_contain_text("0.1")
    expect(bookmark_page.get_by_placeholder("Search for a gene")).to_have_value("")
    expect(
        bookmark_page.get_by_label(
            "View terms in a sortable grid (instead of hierarchical tree)"
        )
    ).to_be_visible()

    bookmark_page.get_by_text("Time series and Gene Selection").click()

    assert_snapshot(
        widget_locator(bookmark_page, "Expression Time Courses")
        .locator("canvas")
        .screenshot(),
        "time_courses_default.png",
    )
    assert_snapshot(
        widget_locator(bookmark_page, "Hierarchical Clustering")
        .locator("canvas")
        .screenshot(),
        "hierarchical_clustering_default.png",
    )
    assert_snapshot(
        widget_locator(bookmark_page, "Differential expressions")
        .locator("canvas")
        .screenshot(),
        "differential_expressions_default.png",
    )


@pytest.mark.parametrize(
    "authenticated", authentication_scenarios(True, False), indirect=True
)
def test_export_visualizations(root_page_with_gene_selection, widgets):
    page: Page = root_page_with_gene_selection

    today = datetime.now()
    prefix = f"{today.strftime('%d_%m_%Y')}_e2e_{today.strftime('%f')}"

    wait_for_widgets_loaded(page, widgets)

    page.get_by_role("button", name="Export").click()
    page.get_by_label("Optional prefix of exported files").fill(prefix)

    with page.expect_download() as download_info:
        page.get_by_role("button", name="Export").click()
        download = download_info.value

        file_name = download.suggested_filename
        assert file_name == f"{prefix}_Report.zip"

        destination_folder = Path("./downloaded_files")
        download.save_as(destination_folder / file_name)

        assert Path(f"downloaded_files/{file_name}").is_file()
