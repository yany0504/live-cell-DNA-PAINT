% Plot traces near synapses
for m_tr=1:length(TraceID)
    figure
    plot3(TraceAll{TraceID(m_tr)}(:,1),TraceAll{TraceID(m_tr)}(:,2),TraceAll{TraceID(m_tr)}(:,3))
end