from pymol import cmd, stored


def _check_color(color):
	assert len(color) == 3
	assert all(0 <= v <= 1 for v in color)
	return tuple(color)


def _blend_colors(a, c1, c2):
	return tuple(v1 * (1 - a) + v2 * a for v1, v2 in zip(c1, c2))


class ColorMap:

	def __init__(self, colors, vmin=0, vmax=1, bad=(0, 0, 0)):

		t_last = -1

		self.colors = []

		for t, color in colors:
			assert t > t_last
			t_last = t
			assert 0 <= t <= 1
			self.colors.append((t, _check_color(color)))

		assert self.colors

		assert vmin < vmax
		self.vmin = vmin
		self.vmax = vmax

		self.bad = _check_color(bad)

	def __call__(self, val):
		val = float(val)

		# Check NaN
		if val != val:
			return self.bad

		# Scale
		val = (val - self.vmin) / (self.vmax - self.vmin)

		min_t, min_color = self.colors[0]
		if val <= min_t:
			return min_color

		last_t = min_t
		last_color = min_color

		for t, color in self.colors[1:]:

			if val < t:
				a = (val - last_t) / (t - last_t)
				return _blend_colors(a, last_color, color)

			last_t = t
			last_color = color

		return color

	def with_range(self, vmin, vmax):
		return ColorMap(self.colors, vmin, vmax, self.bad)


RED_BLUE = ColorMap(
	[
		(0, (1, 0, 0)),
		(.5, (.75, .75, .75)),
		(1, (0, 0, 1)),
	],
)

CMAP = RED_BLUE.with_range(-2, 2)

CMAP2 = ColorMap([
	(0, (.05, .02, .5)),
	(.5, (.8, .28, .47)),
	(1, (.94, .98, .13)),
], bad=(.25, .25, .25))

CMAP3 = ColorMap([
	(0, (.75, .75, .75)),
	#(.5, (1, .75, .75)),
	(1, (1, 0, 0)),
], bad=(.25, .25, .25))


if __name__ == 'pymol':

	cmd.reinitialize()

	# Fetch amyloid structure
	#AMYLOID = cmd.fetch('2N0A', 'amyloid')
	AMYLOID = cmd.fetch('1XQ8', 'amyloid')

	# Display as cartoon with side chains
	cmd.show_as('cartoon', AMYLOID)
	cmd.show('sticks', AMYLOID)

	# Hide all but chain A
	cmd.hide('everything', 'not chain A in bo. ' + AMYLOID)
	cmd.center('chain A in bo. ' + AMYLOID)

	import json

	with open('/Users/student/projects/pubs/structure/pymol_data.json', 'rb') as fobj:
		data = json.loads(fobj.read())

		vals = data['max_charge_mut_delta_fit']

	cmap = CMAP3.with_range(0, max(vals))

	for i, v in enumerate(vals):
		sel = 'resi {} in bo. {}'.format(i + 1, AMYLOID)
		colorname = 'res{}'.format(i + 1)
		cmd.set_color(colorname, list(cmap(v)))
		cmd.color(colorname, sel)

	labels = dict()

	for res, from_aa, to_aa, fit_rel in data['top_charge_muts']:
		label = '{}{}{}: {:+.1f}'.format(from_aa, res + 1, to_aa, fit_rel)
		labels.setdefault(res, []).append(label)

	for res, l in labels.items():
		sel = 'chain A and n. CA and res {} in bo. {}'.format(res + 1, AMYLOID)
		cmd.label(sel, '"' + r'\n'.join(l) + '"')
