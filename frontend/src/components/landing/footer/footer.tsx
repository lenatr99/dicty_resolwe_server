import React, { ReactElement } from 'react';
import {
    Facebook as FacebookIcon,
    Email as EmailIcon,
    Twitter as TwitterIcon,
    LinkedIn as LinkedInIcon,
} from '@mui/icons-material';

import { FooterContainer, FooterWrapper, FooterTitle, SocialNetworkLink } from './footer.styles';
import { SectionContentContainer } from 'components/landing/common/layout.styles';
import SectionNames from 'components/landing/common/constants';

const Footer = (): ReactElement => (
    <FooterContainer>
        <SectionContentContainer id={SectionNames.CONTACT}>
            <FooterWrapper>
                <div>
                    <FooterTitle>Get in touch</FooterTitle>
                    <p>
                        If you have an idea of how to improve dictyExpress
                        <br />
                        or to add extra functionality,{' '}
                        <a href="mailto:info@genialis.com">please contact us.</a>
                    </p>
                </div>
                <div>
                    <div>
                        <SocialNetworkLink href="mailto:info@genialis.com">
                            <EmailIcon />
                        </SocialNetworkLink>
                        <SocialNetworkLink href="https://www.facebook.com/genialisinc/">
                            <FacebookIcon />
                        </SocialNetworkLink>

                        <SocialNetworkLink href="https://twitter.com/Genialis">
                            <TwitterIcon />
                        </SocialNetworkLink>

                        <SocialNetworkLink href="https://www.linkedin.com/company/genialis">
                            <LinkedInIcon />
                        </SocialNetworkLink>
                    </div>
                    <div>
                        <a
                            href="http://www.genialis.com/"
                            rel="noopener noreferrer"
                            target="_blank"
                        >
                            www.genialis.com
                        </a>{' '}
                        {' | '}
                        <a href="mailto:info@genialis.com">info@genialis.com</a>
                    </div>
                    <div>
                        Copyright: Genialis, Inc. {`${new Date().getFullYear()} | `}
                        <a
                            href="https://www.genialis.com/privacy-policy/"
                            rel="noopener noreferrer"
                            target="_blank"
                        >
                            Privacy policy
                        </a>
                        {' | '}
                        <a
                            href="http://www.alexanderwild.com/"
                            rel="noopener noreferrer"
                            target="_blank"
                        >
                            Dictyostelium Image © Alex Wild
                        </a>
                    </div>
                </div>
            </FooterWrapper>
        </SectionContentContainer>
    </FooterContainer>
);

export default Footer;
