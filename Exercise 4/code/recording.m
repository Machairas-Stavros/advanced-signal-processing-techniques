samplingRate = 44100;
bitsPerSample = 16;
numberOfChannels = 1;
recordingDuration = 3;
recObj = audiorecorder(samplingRate,bitsPerSample,numberOfChannels);
recordblocking(recObj,recordingDuration);
audioData = getaudiodata(recObj);
play(recObj)
plot(audioData)