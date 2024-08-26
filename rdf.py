import os, sys
os.environ['OVITO_GUI_MODE'] = '1'
import ovito
import numpy as np
import matplotlib.pyplot as plt

data_file = sys.argv[1]
data_name = os.path.basename(data_file)
data_name_without_extension = os.path.splitext(data_name)[0]

# Load a simulation trajectory consisting of several frames:
pipeline = ovito.io.import_file(data_file)
print("Number of MD frames:", pipeline.source.num_frames)

# Print the list of input particle types.
# They are represented by ParticleType objects attached to the 'Particle Type' particle property.
for t in pipeline.compute().particles.particle_types.types:
    print("Type %i: %s" % (t.id, t.name))

# Insert the RDF calculation modifier into the pipeline:
#pipeline.modifiers.append(ExpressionSelectionModifier(expression = "ParticleType==3 && Position.Z < 100"))
pipeline.modifiers.append(ovito.modifiers.CoordinationAnalysisModifier(cutoff = 10.0, number_of_bins = 100, partial=True))

# Insert the time-averaging modifier into the pipeline, which accumulates
# the instantaneous DataTable produced by the previous modifier and computes a mean histogram.
pipeline.modifiers.append(ovito.modifiers.TimeAveragingModifier(operate_on='table:coordination-rdf'))

# Data export method: Use OVITO's own export function for DataTable objects:
ovito.io.export_file(pipeline, "%s_tot_rdf.txt"%data_name_without_extension, "txt/table", key="coordination-rdf[average]")