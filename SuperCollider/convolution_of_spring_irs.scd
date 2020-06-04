s.boot;

// Uncomment if you want to change the output/input device
//
// ServerOptions.devices;
//
// s.options.outDevice;
// s.options.inDevice;
//
// (
// s.options.outDevice_(
// 	"Built-in Output"
// 	// "Soundflower (2ch)"
// );
// );
//
//
// (
// s.options.inDevice_(
// 	"Built-in Microph"
// 	// "Soundflower (2ch)"
// );
// );




// Change the paths to the correct locations of the impulse response files
(

~fftsize = 2048;

s.waitForBoot {

    {
        var ir, irbuffer, bufsize;
        // // MONO ONLY
        // pre-existing impulse response sound files

		// irbuffer = Buffer.read(s,"/Users/mariusonofrei/Documents/Physical-Modelling-of-Spring-Reverberation-master/Helical_Spring_IR_Database/irEstimate_K_0p02_q_800_gamma_1800.wav");
		irbuffer = Buffer.read(s,"/Users/mariusonofrei/Documents/Physical-Modelling-of-Spring-Reverberation-master/Helical_Spring_IR_Database/irEstimate_K_0p16_q_800_gamma_1800.wav");
		s.sync;

        bufsize = PartConv.calcBufSize(~fftsize, irbuffer);

        ~irspectrum = Buffer.alloc(s, bufsize, 1);
        ~irspectrum.preparePartConv(irbuffer, ~fftsize);
        ~irspectrum = ~irspectrum.normalize;

        s.sync;

        irbuffer.free; // don't need time domain data anymore, just needed spectral version
    }.fork;

}
)


(
SynthDef.new(\soundin_synth,{
	arg gate=0,noiseHz=0.1,amp=1,bufnum=~irspectrum.bufnum,freq=440,rq=1,bfmix=0,pan=0;
	var sig;
	sig = SoundIn.ar(0);
	sig = PartConv.ar(sig, ~fftsize, bufnum);
	// sig = LeakDC.ar(sig, 0.99);
	sig = XFade2.ar(sig,BPF.ar(sig,freq,rq,1/rq.sqrt),bfmix*2-1);
	sig = Pan2.ar(sig!2,pan,amp);
	sig = Limiter.ar(sig,0.4,0.01);
	// sig = sig.tanh;
	Out.ar(0,sig);
}).add;
)


s.scope;

(
x = Synth.new(\soundin_synth,[
	\bufnum,~irspectrum.bufnum,
	\gate,1,
	\amp,1,
	\bfmix,0.0,
	\rq,1/150,
	\freq,1000]);
)

x.set(\bfmix,0.99,\freq,110)

x.free;








