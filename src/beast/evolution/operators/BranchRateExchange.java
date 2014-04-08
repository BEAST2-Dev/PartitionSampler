package beast.evolution.operators;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.QuantileUCRelaxedClock;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import org.apache.commons.math.MathException;

/**
 * @author Chieh-Hsi Wu
 */
public class BranchRateExchange extends Exchange{
    public Input<RealParameter> quantilesInput = new Input<RealParameter>(
            "quantiles",
            "the quantiles of the prior distribution on the rates associated with nodes in the tree for sampling of individual rates among branches.",
            Input.Validate.REQUIRED
    );
    public Input<QuantileUCRelaxedClock> quantileUCRelaxedClockInput = new Input<QuantileUCRelaxedClock>(
            "quantileUCRelaxedClock",
            "The relaxed clock used to model temporal heterogeneity in rates.",
            Input.Validate.REQUIRED
    );

    private RealParameter quantiles;
    private QuantileUCRelaxedClock quantileUCRelaxedClock;
    public void initAndValidate(){
        quantiles = quantilesInput.get();
        quantileUCRelaxedClock = quantileUCRelaxedClockInput.get();

        super.initAndValidate();
    }

    boolean goodQuantileProposal = true;

    public double narrow(final Tree tree) {

        double logq = super.narrow(tree);
        if(goodQuantileProposal){
            return logq;
        }else{
            return Double.NEGATIVE_INFINITY;
        }
    }

    public double wide(final Tree tree) {

        double logq = super.wide(tree);
        if(goodQuantileProposal){
            return logq;
        }else{
            goodQuantileProposal = true;
            return Double.NEGATIVE_INFINITY;
        }
    }


    protected void exchangeNodes(Node i, Node j,
                                 Node iP, Node jP) {
        try{
            double oldiSubst = i.getLength()*quantileUCRelaxedClock.getRawRateForBranch(i);
            double oldjSubst = j.getLength()*quantileUCRelaxedClock.getRawRateForBranch(j);

            // precondition iP -> i & jP -> j
            replace(iP, i, j);
            replace(jP, j, i);
            // postcondition iP -> j & iP -> i

            double newiQ = quantileUCRelaxedClock.getQuantile(oldiSubst/i.getLength());
            double newjQ = quantileUCRelaxedClock.getQuantile(oldjSubst/j.getLength());
            if(newiQ == 1.0 || newjQ == 1.0){
                throw new MathException();
            }
            quantiles.setValue(quantileUCRelaxedClock.getQuantileIndex(i),newiQ);
            quantiles.setValue(quantileUCRelaxedClock.getQuantileIndex(j),newjQ);
        }catch(MathException e){
            goodQuantileProposal = false;

        }
    }
}
