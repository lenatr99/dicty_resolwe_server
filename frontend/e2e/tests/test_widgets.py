import pytest
from playwright.sync_api import Locator, Page, expect
from tests.utils import (
    add_genes_to_selection,
    authentication_scenarios,
    wait_for_widgets_loaded,
    widget_locator,
)


def assert_tree_items(widget: Locator, tree_items: list[str]):
    for item in tree_items[:-1]:
        expect(widget.get_by_role("gridcell", name=item, exact=True)).to_be_visible()
        expect(
            widget.get_by_role("gridcell", name=item, exact=True).get_by_test_id(
                "ArrowDropDownIcon"
            )
        ).to_be_visible()
    expect(
        widget.get_by_role("gridcell", name=tree_items[-1], exact=True)
    ).to_be_visible()


@pytest.mark.parametrize(
    "authenticated", authentication_scenarios(True, False), indirect=True
)
def test_gene_selection(root_page_with_gene_selection, genes, widgets, assert_snapshot):
    page: Page = root_page_with_gene_selection

    wait_for_widgets_loaded(page, widgets)

    assert_snapshot(
        page.get_by_test_id("genes-expressions-line-chart")
        .locator("canvas")
        .screenshot(),
        "time_courses_default.png",
    )

    page.get_by_role("button", name="enoB").click()
    expect(page.get_by_role("heading", name="Gene information")).to_be_visible()
    expect(page.get_by_text("DDB_G0268214")).to_be_visible()
    expect(page.get_by_text("enolase B")).to_be_visible()

    with page.expect_popup() as dictybase_popup_info:
        page.get_by_role("link", name="Open in dictyBaseOpen in dictyBase.").click()
    dictybase_popup_page = dictybase_popup_info.value
    expect(dictybase_popup_page).to_have_url("http://dictybase.org/gene/DDB_G0268214")
    dictybase_popup_page.close()

    page.get_by_role("button", name="Highlight").click()
    page.get_by_text("Expression Time Courses").click()

    assert_snapshot(
        page.get_by_test_id("genes-expressions-line-chart")
        .locator("canvas")
        .screenshot(),
        "time_courses_highlighted.png",
    )

    assert_snapshot(
        widget_locator(page, "Hierarchical Clustering").locator("canvas").screenshot(),
        "hierarchical_clustering_highlighted.png",
    )

    page.get_by_test_id("genes-expressions-line-chart").click()
    page.get_by_text("Expression Time Courses").click()

    assert_snapshot(
        page.get_by_test_id("genes-expressions-line-chart")
        .locator("canvas")
        .screenshot(),
        "time_courses_default.png",
    )

    page.get_by_label("Copy 4 genes to clipboard").click()
    clipboard_content = page.evaluate("navigator.clipboard.readText()")
    assert clipboard_content == ", ".join(genes)

    page.get_by_role("button", name=genes[0]).get_by_test_id("CancelIcon").click()
    expect(page.get_by_role("button", name=genes[0])).not_to_be_visible()


@pytest.mark.parametrize(
    "authenticated", authentication_scenarios(True, False), indirect=True
)
def test_gene_history(root_page, genes):
    page: Page = root_page

    history_genes = ["pks6", "pkaR", "pks5", "pkgB", "ppk1", "pkaC"]
    add_genes_to_selection(page, history_genes)

    page.get_by_label("Save current gene set to local storage").click()
    expect(page.get_by_text("Gene set saved.")).to_be_visible()

    page.get_by_label("Clear all").click()
    for gene in history_genes:
        expect(page.get_by_role("button", name=gene)).not_to_be_visible()

    add_genes_to_selection(page, genes)

    page.get_by_role("button", name="History").click()
    expect(page.get_by_role("heading", name="Gene List History")).to_be_visible()

    page.get_by_text(", ".join(history_genes)).click()
    for gene in history_genes:
        expect(page.get_by_role("button", name=gene)).to_be_visible()


@pytest.mark.parametrize(
    "authenticated", authentication_scenarios(True, False), indirect=True
)
def test_expression_time_courses(root_page_with_gene_selection, assert_snapshot):
    page: Page = root_page_with_gene_selection

    wait_for_widgets_loaded(page, ["Expression Time Courses"])

    widget: Locator = widget_locator(page, "Expression Time Courses")

    assert_snapshot(
        widget.locator("canvas").screenshot(),
        "time_courses_default.png",
    )

    widget.get_by_role("button", name="Find similar genes").click()

    expect(page.get_by_role("heading", name="Find Similar Genes")).to_be_visible()

    expect(page.get_by_label("Gene", exact=True)).to_contain_text("aplA")
    expect(page.get_by_text("No Rows To Show")).to_be_visible()

    page.get_by_role("button", name="Find").click()

    expect(page.get_by_role("progressbar")).not_to_be_visible()
    expect(page.get_by_text("Loading...")).not_to_be_visible()
    expect(page.get_by_text("No Rows To Show")).not_to_be_visible()

    page.get_by_role("textbox").fill("ctps")
    expect(page.get_by_role("gridcell", name="CTP synthase")).to_be_visible()

    page.get_by_text("63.59726ctpsCTP synthaseDDB_G0280567").get_by_role(
        "checkbox"
    ).check()

    page.get_by_role("button", name="Select", exact=True).click()

    assert_snapshot(
        widget.locator("canvas").screenshot(),
        "time_courses_similar_genes_default.png",
    )

    widget.get_by_label("Legend").check()

    assert_snapshot(
        widget.locator("canvas").screenshot(),
        "time_courses_similar_genes_legend.png",
    )

    widget.get_by_role("button", name="Compare to (0 experiments)").click()
    expect(
        page.get_by_role("heading", name="Select Time series To Compare")
    ).to_be_visible()

    page.get_by_role("textbox").fill("transcriptome")
    expect(
        page.get_by_role("gridcell", name="05. lncRNA transcriptome")
    ).to_be_visible()

    page.get_by_text(
        "05. lncRNA transcriptomeD. discoideumAX4Filter DevelopmentHL5Rosengarten et. al."
    ).get_by_role("checkbox").check()

    page.get_by_role("button", name="Close").click()

    assert_snapshot(
        widget.locator("canvas").screenshot(),
        "time_courses_compare_series_legend.png",
    )

    widget.get_by_label("Color by time series").check()

    assert_snapshot(
        widget.locator("canvas").screenshot(),
        "time_courses_compare_series_color_by_time_series_legend.png",
    )

    widget.get_by_label("Legend").uncheck()

    assert_snapshot(
        widget.locator("canvas").screenshot(),
        "time_courses_compare_series_color_by_time_series_no_legend.png",
    )


@pytest.mark.parametrize(
    "authenticated", authentication_scenarios(True, False), indirect=True
)
def test_gene_oncology_enrichment(root_page_with_gene_selection, genes):
    page: Page = root_page_with_gene_selection

    wait_for_widgets_loaded(page, ["Gene Ontology Enrichment"])

    widget: Locator = widget_locator(page, "Gene Ontology Enrichment")

    expect(
        page.get_by_text("Computing Gene Ontology Enrichment for 4 genes.")
    ).not_to_be_visible()

    tree_items = [
        "biological adhesion",
        "adhesion of symbiont to host",
        "adhesion of symbiont to host cell",
        "virion attachment to host cell",
        "receptor-mediated virion attachment to host cell",
        "receptor binding",
        "G-protein coupled receptor binding",
    ]

    assert_tree_items(widget, tree_items)

    widget.get_by_role("gridcell", name=tree_items[0]).get_by_test_id(
        "ArrowDropDownIcon"
    ).click()

    for item in tree_items[1:]:
        expect(
            widget.get_by_role("gridcell", name=item, exact=True)
        ).not_to_be_visible()

    widget.get_by_role("gridcell", name=tree_items[0]).get_by_test_id(
        "ArrowDropDownIcon"
    ).click()

    widget.get_by_label(
        "View terms in a sortable grid (instead of hierarchical tree)"
    ).click()

    for item in tree_items:
        expect(widget.get_by_role("gridcell", name=item, exact=True)).to_be_visible()
        expect(
            widget.get_by_role("gridcell", name=item, exact=True).get_by_test_id(
                "ArrowDropDownIcon"
            )
        ).not_to_be_visible()

    widget.get_by_label("View terms in a hierarchical tree").click()

    assert_tree_items(widget, tree_items)

    expect(widget.get_by_label("Aspect")).to_contain_text("Biological process")
    widget.get_by_label("Aspect").click()

    expect(page.get_by_role("option", name="Cellular component")).to_be_visible()
    expect(page.get_by_role("option", name="Molecular function")).to_be_visible()

    page.get_by_role("option", name="Cellular component").click()

    expect(widget.get_by_label("p-value")).to_contain_text("0.1")

    widget.get_by_label("p-value").click()
    page.get_by_role("option", name="0.05").click()

    expect(
        widget.get_by_text("Computing Gene Ontology Enrichment for 4 genes.")
    ).to_be_visible(timeout=10000)

    wait_for_widgets_loaded(widget, ["Gene Ontology Enrichment"])

    expect(widget.get_by_text("Enriched terms not found")).not_to_be_visible()


@pytest.mark.parametrize(
    "authenticated", authentication_scenarios(True, False), indirect=True
)
def test_hierarchical_clustering(
    root_page_with_gene_selection, widgets, assert_snapshot
):
    page: Page = root_page_with_gene_selection

    wait_for_widgets_loaded(page, widgets)

    widget: Locator = widget_locator(page, "Hierarchical Clustering")

    assert_snapshot(
        widget.locator("canvas").screenshot(),
        "hierarchical_clustering_default.png",
    )

    expect(widget.get_by_label("Clustering Linkage")).to_contain_text("Average")

    widget.get_by_label("Clustering Linkage").click()
    expect(page.get_by_role("option", name="Single")).to_be_visible()
    expect(page.get_by_role("option", name="Average")).to_be_visible()
    expect(page.get_by_role("option", name="Complete")).to_be_visible()

    page.get_by_role("option", name="Complete").click()
    expect(widget.get_by_label("Clustering Linkage")).to_contain_text("Complete")

    expect(widget.get_by_test_id("ProgressBar")).not_to_be_visible()

    page.wait_for_timeout(1000)

    assert_snapshot(
        widget.locator("canvas").screenshot(),
        "hierarchical_clustering_complete.png",
    )

    expect(widget.get_by_label("Distance Measure")).to_contain_text("Euclidean")

    widget.get_by_label("Distance Measure").click()
    expect(page.get_by_role("option", name="Euclidean")).to_be_visible()
    expect(page.get_by_role("option", name="Pearson")).to_be_visible()
    expect(page.get_by_role("option", name="Spearman")).to_be_visible()

    page.get_by_role("option", name="Pearson").click()
    expect(widget.get_by_label("Distance Measure")).to_contain_text("Pearson")

    expect(widget.get_by_test_id("ProgressBar")).not_to_be_visible()

    page.wait_for_timeout(1000)

    assert_snapshot(
        widget.locator("canvas").screenshot(),
        "hierarchical_clustering_complete_pearson.png",
    )


@pytest.mark.parametrize(
    "authenticated", authentication_scenarios(True, False), indirect=True
)
def test_differential_expressions(root_page_with_gene_selection, assert_snapshot):
    page: Page = root_page_with_gene_selection

    wait_for_widgets_loaded(page, ["Differential expressions"])

    widget: Locator = widget_locator(page, "Differential expressions")
    differential_expressions_plot: Locator = widget.locator("canvas")

    assert_snapshot(
        differential_expressions_plot.screenshot(),
        "differential_expressions_default.png",
    )

    differential_expressions_plot.hover(position={"x": 330, "y": 230})
    page.mouse.down()
    differential_expressions_plot.hover(position={"x": 250, "y": 150})
    page.mouse.up()

    expect(
        page.get_by_role("heading", name="Selected Differential Expression Genes")
    ).to_be_visible()

    page.get_by_role("textbox").fill("ugt2")
    expect(page.get_by_role("gridcell", name="ugt2")).to_be_visible()

    expect(page.get_by_label("Append selected genes to Genes module")).to_be_checked()

    row = page.get_by_text("DDB_G0268540ugt23.6110002342.4055757351295775")
    row.get_by_role("checkbox").check()
    page.get_by_role("button", name="Select", exact=True).click()

    assert_snapshot(
        differential_expressions_plot.screenshot(),
        "differential_expressions_added_gene.png",
    )
