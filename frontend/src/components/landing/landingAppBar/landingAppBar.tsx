import React, { ReactElement } from 'react';
import { useNavigate } from 'react-router-dom';
import { MenuItem } from '@mui/material';
import { ShortDescription, LandingAppBarWrapper, NavButton } from './landingAppBar.styles';
import SectionNames from 'components/landing/common/constants';
import { scrollToTargetAdjusted } from 'utils/documentHelpers';
import { appBarHeight } from 'components/app/globalStyle';
import { ContentContainer } from 'components/landing/common/layout.styles';
import DictyAppBar from 'components/common/dictyAppBar/dictyAppBar';

const scrollButtonsData = [
    {
        name: 'Home',
        sectionHash: '',
    },
    {
        name: 'Features',
        sectionHash: SectionNames.FEATURES,
    },
    {
        name: 'Screenshots',
        sectionHash: SectionNames.SCREENSHOTS,
    },
    {
        name: 'Citing',
        sectionHash: SectionNames.CITING,
    },
    {
        name: 'Contact',
        sectionHash: SectionNames.CONTACT,
    },
];

const LandingAppBar = (): ReactElement => {
    const navigate = useNavigate();

    const navigateToDictyExpress = (): void => {
        navigate('/bcm');
    };

    const scrollToSection = (hash: string): void => {
        scrollToTargetAdjusted(hash, appBarHeight);
    };

    const getScrollButton = (name: string, sectionHash: string): ReactElement => (
        <NavButton key={name} onClick={(): void => scrollToSection(sectionHash)}>
            {name}
        </NavButton>
    );

    const runDictyButton = <NavButton onClick={navigateToDictyExpress}>Run dictyExpress</NavButton>;

    const desktopSection = (
        <>
            {scrollButtonsData.map((buttonData) =>
                getScrollButton(buttonData.name, buttonData.sectionHash),
            )}
            {runDictyButton}
            <ShortDescription>
                Analysis of gene expression in <i>Dictyostelium</i>
            </ShortDescription>
        </>
    );

    const mobileMenu = (
        <div>
            {scrollButtonsData.map((buttonData) => (
                <MenuItem key={buttonData.name}>
                    {getScrollButton(buttonData.name, buttonData.sectionHash)}
                </MenuItem>
            ))}
            <MenuItem>{runDictyButton}</MenuItem>
        </div>
    );

    return (
        <>
            <LandingAppBarWrapper position="sticky">
                <ContentContainer>
                    <DictyAppBar desktopSection={desktopSection} mobileMenu={mobileMenu} />
                </ContentContainer>
            </LandingAppBarWrapper>
        </>
    );
};

export default LandingAppBar;
