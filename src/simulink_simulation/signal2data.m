function data = signal2data(dataset, signal)
    signal = getElement(dataset, signal).Values;

    data.Data = squeeze(signal.Data);
    data.Time = signal.Time;

end