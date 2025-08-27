# E2e test automation

This folder holds scripts for testing the dictyExpress frontend.
All scripts should only be run with **basic account**!

## Prerequisites

Prepare a basic user on Expressions. The user should have staff and superuser status set to `false`, and it should be assigned the permission set `e2e_test`.

Before running the tests, you should set the following environment variables with your prepared user login info:

-   `EXPRESSIONS_E2E_USERNAME` (basic user which is used for all tests)
-   `EXPRESSIONS_E2E_PASSWORD`

You should also install all required libraries:

```bash
pip install -r requirements.txt
```

Install the required browsers for Playwright:

```bash
python -m playwright install
```

## Running tests

Tests can be ran in different environments:

-   `"app"` runs the tests on app.dictyexpress.org
-   `"qa"` runs the tests on qa.dictyexpress.org
-   `"dev"` runs the tests on localhost

To run the tests you must specify the `--env` argument.

Inside the folder `e2e` run:

```bash
python -m pytest tests --env="qa"
```

To skip public user tests run the command with `"authenticated"` mark:

```bash
python -m pytest tests --env="app" -m authenticated
```

To run a single test or test file provide the full relative path to the test or test file.
For example, to only run the `test_page_layout` test:

```bash
python -m pytest tests/general/test_general.py::test_page_layout --env="qa"
```

To see the verbose test logs during execution, run:

```bash
python -m pytest tests --env="qa" -v -s
```

There are two different modes in which you can run tests, **headed** and **headless** mode.
In headed mode the tests open browser windows, while in headless the browser is opened only in the background.
The default mode is headless. To see what is happening during tests and increase the timeout between each step, run:

```bash
python -m pytest tests --env="qa" --slowmo=1000 --headed
```

To [run tests in debug mode](https://playwright.dev/python/docs/running-tests#debugging-tests), set the environment variable `PWDEBUG` to `1`.

You can also [record a trace] of a test. To do so set the environment variable `CI` to `true`. Traces will be recorded for each test ran in the `e2e/artifacts/traces` folder.
Then to view the trace, run the following command:

```bash
python -m playwright show-trace "artifacts/traces/[trace_file_name]"
```

If you want to create html report of the tests, add
parameters `--html=artifacts/reports/report.html` and `--self-contained-html`
to the command above.

To generate a verification checklist for Confluence, use the `verification_checklist.py` script.
This script uses the `report.html` file as input and outputs a DOCX file in the same folder as the HTML file:

```bash
python verification_checklist.py --file-path "artifacts/reports/report.html" --output-fn "verification_checklist.docx"
```
