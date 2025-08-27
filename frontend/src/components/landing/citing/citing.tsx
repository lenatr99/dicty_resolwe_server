import React, { ReactElement } from 'react';
import { SectionContentContainer } from 'components/landing/common/layout.styles';
import SectionNames from 'components/landing/common/constants';
import { Title } from 'components/landing/common/title.styles';

const Citing = (): ReactElement => (
    <SectionContentContainer id={SectionNames.CITING} $centerText $paddingBottom={35}>
        <Title>Citing dictyExpress</Title>
        <p>
            Stajdohar M, Rosengarten RD, Kokosar J, Jeran L, Blenkus D, Shaulsky G, Zupan B.{' '}
            <a href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1706-9">
                &quot;dictyExpress: a web-based platform for sequence data management and analytics
                in <i>Dictyostelium</i> and beyond.&quot;
            </a>{' '}
            BMC Bioinformatics. 2017. 18(1):291
        </p>
    </SectionContentContainer>
);

export default Citing;
