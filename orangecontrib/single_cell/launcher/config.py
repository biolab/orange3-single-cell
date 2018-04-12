from types import SimpleNamespace

import pkg_resources

FEEDBACK_URL = "http://orange.biolab.si/survey/long.html"


class WelcomeScreenSpecs(SimpleNamespace):
    class Item(SimpleNamespace):
        path = ""  # type: str
        icon = ""  # type: str
        active_icon = ""  # type: str
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
    def workflow_filename(filename):
        return pkg_resources.resource_filename(
            "orangecontrib.single_cell.tutorials", filename)

    def image_filename(filename):
        return pkg_resources.resource_filename(
            "orangecontrib.single_cell.launcher", "icons/"+filename)

    items = [
        WelcomeScreenSpecs.Item(
            path=workflow_filename("Showcase-CellType.ows"),
            icon=image_filename("CellTypeDiscovery.svg"),
            active_icon=image_filename("CellTypeDiscoveryColor.svg"),
            text="Cell Type Discovery",
            tooltip="Use Louvain clustering and tSNE to find cell types",
        ),
        WelcomeScreenSpecs.Item(
            path=workflow_filename("Showcase-Markers.ows"),
            icon=image_filename("Biomarkers.svg"),
            active_icon=image_filename("BiomarkersColor.svg"),
            text="Biomarkers",
            tooltip="Visualize markers on tSNE plot"
        ),
        WelcomeScreenSpecs.Item(
            path=workflow_filename("Showcase-Prediction.ows"),
            icon=image_filename("PopulationPrediction.svg"),
            active_icon=image_filename("PopulationPredictionColor.svg"),
            text="Population Prediction",
            tooltip="Predict population with some widgets"
        )
    ]
    return WelcomeScreenSpecs(
        image=image_filename("welcome-background.png"),
        css="""
            StartItem {
                background-color: rgb(214, 232, 236);
                color: rgb(109, 172, 178);
            }

            StartItem::item {
                border: 1px solid transparent;
                padding-top: 1.5em;
            }
            """,
        items=items,
    )
