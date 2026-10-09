function phasered = phaseReduction(phase)
    
    phasered = mod(phase,2*pi);

    pos = find(phase>=0);
    neg = find(phase<0);
    pos_pito2pi = find(phasered(pos)>pi);
    neg_pito2pi = find(phasered(neg)>pi);

    phasered(pos(pos_pito2pi)) = phasered(pos(pos_pito2pi))-2*pi;
    phasered(neg(neg_pito2pi)) = phasered(neg(neg_pito2pi))-2*pi;

end