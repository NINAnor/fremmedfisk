import os
project_instance = QgsProject.instance()
manager = project_instance.layoutManager()
outdir = '/home/stefan/fremmedfisk/plots'

for layout in manager.layouts():
    print("Exporting {}...".format(layout.name()))
    exporter = QgsLayoutExporter(layout)
    exporter.ImageExportSettings().cropToContents = True
    exporter.ImageExportSettings().dpi = 300
    outpath = os.path.join(outdir, layout.name() + ".png")
    exporter.exportToImage(outpath, exporter.ImageExportSettings())
