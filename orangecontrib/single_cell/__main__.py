from Orange.canvas.__main__ import main
from Orange.canvas import registry as orange_registry


wd = orange_registry.WidgetDiscovery.widget_description
def widget_description(self, module, widget_name=None,
                       category_name=None, distribution=None):
    """Move tSNE widget to Visualize category

    A better way to do this would be to change the behaviour of the
    widget_description method to *not* overwrite the category when
    manually specified on widget.
    """

    desc = wd(
        self, module, widget_name, category_name, distribution)

    if desc.qualified_name == 'orangecontrib.single_cell.widgets.owtsne.OWtSNE':
        desc.category = 'Visualize'
    return desc

orange_registry.WidgetDiscovery.widget_description = widget_description


if __name__ == "__main__":
    main()
