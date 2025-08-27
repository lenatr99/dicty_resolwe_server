import React, { ReactElement } from 'react';
import { CreditsLogos } from './credits.styles';
import SectionNames from 'components/landing/common/constants';
import {
    SectionContentContainer,
    DarkSectionContainer,
} from 'components/landing/common/layout.styles';
import { Title } from 'components/landing/common/title.styles';
import baylorLogo from 'images/bcmLogo.png';
import friLogo from 'images/ulfriLogo.png';
import genialisLogo from 'images/genialisLogo.png';

const Credits = (): ReactElement => (
    <DarkSectionContainer>
        <SectionContentContainer id={SectionNames.CREDITS} $centerText>
            <Title>Credits</Title>
            <p>
                dictyExpress was built in collaboration with
                <br />
                <a href="http://www.biolab.si">
                    Bioinformatics Laboratory at University of Ljubljana
                </a>
                <br />
                <a href="https://www.bcm.edu/people/view/b17c52a8-ffed-11e2-be68-080027880ca6">
                    Gad Shaulsky&apos;s
                </a>{' '}
                and{' '}
                <a href="https://www.bcm.edu/people/view/b1571132-ffed-11e2-be68-080027880ca6">
                    Adam Kuspa&apos;s
                </a>{' '}
                labs at Baylor College of Medicine
            </p>
            <p>
                Powered by <a href="http://www.genialis.com">Genialis</a>
            </p>
            <CreditsLogos>
                <a href="https://www.bcm.edu/">
                    <img src={baylorLogo} alt="Baylor College of Medicine" />
                </a>
                <a href="http://www.fri.uni-lj.si/">
                    <img src={friLogo} alt="University of Ljubljana" />
                </a>
                <a href="http://www.genialis.com/">
                    <img src={genialisLogo} alt="Genialis" />
                </a>
            </CreditsLogos>
        </SectionContentContainer>
    </DarkSectionContainer>
);

export default Credits;
