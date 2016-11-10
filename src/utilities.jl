"""
  _downshift(G,L,shift1,shift2,i,N,cut)

Shift the elements of `L` into `G` at indices `i:N`. Elements below index `cut`
are shifted `shift1` steps and above `cut` are shifted `shift2` steps.
"""
function _downshift!(G,L,shift1,shift2,i,N,cut)
  if i <= cut  # downshift first part of G
    _downshift!(view(G,i:cut), view(L,i:cut), shift1)
  end          # downshift second part of G
  _downshift!(view(G,max(cut+1,i):N), view(L,max(cut+1,i):N), shift2)
end

function _downshift!(G,L,shift)
  if length(G) < shift
    for i in eachindex(G)
      G[i] = zero(eltype(G))
    end
    return
  end
  G[shift+1:end] = L[1:end-shift]
  G[1:shift] = zeros(shift)
end
