from types import SimpleNamespace

import pkg_resources

FEEDBACK_URL = "http://orange.biolab.si/survey/long.html"


class WelcomeScreenSpecs(SimpleNamespace):
    class Item(SimpleNamespace):
        path = ""  # type: str
        icon = ""  # type: str
        text = ""  # type: str
        tooltip = ""  # type: str

    image = ""  # type: str
    css = ""  # Any extra css to apply to the welcome screen widget
    items = []  # type: List[WelcomeScreenSpecs.Item]


def welcome_screen_specs():
    """
    Returns
    -------
    spec : WelcomeScreenSpecs
    """
    def resource_filename(filename):
        return pkg_resources.resource_filename(
            "Orange.canvas.application.workflows", filename)

    def image_filename(filename):
        return pkg_resources.resource_filename(
            "orangecontrib.single_cell.launcher", "icons/"+filename)

    items = [
        WelcomeScreenSpecs.Item(
            path=resource_filename("110-file-and-data-table-widget.ows"),
            icon=image_filename("CellTypeDiscovery.svg"),
            text="Cell Type Discovery",
            tooltip="Something to do with latin",
        ),
        WelcomeScreenSpecs.Item(
            path=resource_filename("310-clustering.ows"),
            icon=image_filename("Biomarkers.svg"),
            text="Biomarkers",
            tooltip="And now for something completely different"
        ),
        WelcomeScreenSpecs.Item(
            path=resource_filename("450-cross-validation.ows"),
            icon=image_filename("PopulationPrediction.svg"),
            text="Population Prediction",
            tooltip=(
                '<p style="color: #555;">'
                '<div style="font-size: large; white-space: nowrap;" >'
                'And now for something <b>completely different</b><hr/>'
                '</div>'
                '&hellip;'
                '</p>'
            )
        )
    ]
    return WelcomeScreenSpecs(
        image=image_filename("welcome-background.png"),
        css="""
            StartItem { background-color: rgb(214, 232, 236) }
            """,
        items=items,
    )
