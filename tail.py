class Tail(object):
    def __init__(self, velocity, temperature):
        self.temperature = temperature
        self.velocity = velocity
        kine_vis = 1.34 * 10 ** - 5. + 9.31477 * 10 ** - 8. * self.temperature
        self.airDensity = 1.28912 - 0.004122391 * self.temperature
        self.dynpres = 0.5 * self.airDensity * self.velocity ** 2.0

    def calc_htaildrag(self, htailArea):
        CDhtail = 0.008
        self.htailDrag = CDhtail * htailArea * self.dynpres
        self.htailPower = self.htailDrag * self.velocity

    def calc_vtaildrag(self, vtailArea):
        CDvtail = 0.008
        self.vtailDrag = CDvtail * vtailArea * self.dynpres
        self.vtailPower = self.vtailDrag * self.velocity