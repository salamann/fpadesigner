class Body(object):
    def __init__(self, velocity, temperature):
        self.temperature = temperature
        self.velocity = velocity
        kine_vis = 1.34 * 10 ** - 5. + 9.31477 * 10 ** - 8. * self.temperature
        self.airDensity = 1.28912 - 0.004122391 * self.temperature
        self.dynpres = 0.5 * self.airDensity * self.velocity ** 2.0

    def fairdragcalc(self, fairArea):
        """CD=0.10 : F-TEC technical report"""
        """CD=0.20 : T-MIT windtunnel data"""
        CDfair = 0.15
        self.fairringDrag = CDfair * fairArea * self.dynpres
        self.fairringPower = self.fairringDrag * self.velocity

    def framedragcalc(self, framearea):
        CDframe = 0.1
        self.frameDrag = CDframe * framearea * self.dynpres
        self.framePower = self.frameDrag * self.velocity


if __name__ == '__main__':
    pass
