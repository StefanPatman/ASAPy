from itaxotools.common.param import Field, Group

def params():
	return Group(key='root', children=[
		Group(key='general', label='General', children=[
			Field(key='all',
				label='Generate all files',
				doc=("Generate all partition files."),
				type=bool,
				default=True
			),
			Field(
				key = 'mega',
				label = 'MEGA CSV',
				doc = ("Set if the provided distance Matrix is CSV file\nexported from MEGA 6 or X."),
				type = bool,
				default = False
			),
			Field(
				key = 'sequence_length',
				label = 'Sequence length',
				doc = ("Original length of the sequence."),
				type = int,
				default = 600
			),
		]),
		Group(key='advanced', label='Advanced', children=[
			Field(key='number',
				label='Scores Kept',
				doc=("Number of results with the highest scores to be displayed."),
				type=int,
				default=10
			),
			Field(key='seuil_pvalue',
				label='Probability',
				doc=("Limit for results to be reported."),
				type=float,
				default=0.01
			),
			Field(key='seed',
				label='Seed',
				doc=("Use fixed seed value. If you don’t want to use a fixed seed value, set to -1."),
				type=int,
				default=-1
			),
		]),
		Group(key='distance', label='Distance', children=[
			Field(key='method',
				label='Method',
				doc=("If you provide a fasta file, you can select the substitution model to compute the distances."),
				type=int,
				list={
                    0: 'Kimura-2P (K80)',
                    1: 'Jukes-Cantor (JC69)',
                    2: 'Tamura-Nei (TN93)',
                    3: 'Simple Distance',
				},
				default=1
			),
			Field(key='rate',
				label='Kimura TS/TV',
				doc=("Transition/transversion for Kimura 3-P distance."),
				type=float,
				default=2.0
			),
		]),
	])
